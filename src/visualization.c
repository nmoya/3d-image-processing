#include "visualization.h"

/* ------------------------------- Private Methods -------------------------*/

/* Methods used for surface and volume rendering */
float     PhongShading(GraphicalContext *gc, int p, Vector N, float dist);
float     DepthShading(GraphicalContext *gc, int p, Vector N, float dist);
Vector    ObjectNormalOnTheFly(GraphicalContext *gc, Image *image, int po, AdjRel *A);
void      RenderFromSRBuffers(GraphicalContext *gc, Image *image);
void      ResetSRBuffers(GraphicalContext *gc);
float     OpacityValueAtPoint(GraphicalContext *gc, Point P);
Vector    NormalVectorAtPoint(GraphicalContext *gc, Point P);
void      SurfaceShadingAlongRay(GraphicalContext *gc, Point P0, Point P1, Point Pn, float *red, float *green, float *blue);
void      VoxelInfoByRayCasting(GraphicalContext *gc, Point P0, Point P1, Point Pn, int *voxel, float *depth, int *object);
void      VolumeShadingAlongRay(GraphicalContext *gc, Point P0, Point P1, Point Pn, float *value);
void      VoxelSplatting(GraphicalContext *gc, Image *image, Point P, int p, float phong_val, float opac);
Image     *SurfaceRenderingByRayCasting(GraphicalContext *gc);
Image     *SurfaceRenderingBySplatting(GraphicalContext *gc);
Image     *VolumeRenderingByRayCasting(GraphicalContext *gc);

/* Methods to create the graphical context */

void              SetFTBaccess(GraphicalContext *gc);
PhongModel       *CreatePhongModel(FImage *scene);
SRBuffers        *CreateSRBuffers(Image *label);
ObjectAttributes *CreateObjectAttributes(Image *label, int *nobjs);
int               GetNormalIndex(Vector normal);
Vector           *CreateNormalTable();
void              SetSceneFaces(GraphicalContext *gc);


float PhongShading(GraphicalContext *gc, int p, Vector N, float dist)
{
    float cos_theta;
    float cos_2theta, pow, phong_val = 0.0;

    cos_theta  = VectorInnerProduct(gc->viewdir->V, N);
    if (cos_theta > 1.0)
        cos_theta = 1.0;

    if (cos_theta > Epsilon)   /* |angle| <= 90° */
    {

        cos_2theta = 2 * cos_theta * cos_theta - 1;

        if (cos_2theta <= Epsilon)   /* |angle| >= 45° */
        {
            pow = 0.;
        }
        else
        {
            pow = 1.;
            for (int k = 0; k < gc->phong->ns; k++)
                pow = pow * cos_2theta;
        }

        phong_val = gc->phong->ka + gc->phong->Idist[(int)dist] * (gc->phong->kd * cos_theta + gc->phong->ks * pow);

    }

    return (phong_val);
}

float DepthShading(GraphicalContext *gc, int p, Vector N, float dist)
{
    float cos_theta;
    float depth_val = 0.0;

    cos_theta  = VectorInnerProduct(gc->viewdir->V, N);
    if (cos_theta > 1.0)
        cos_theta = 1.0;

    if (cos_theta > Epsilon)   /* |angle| <= 90° */
    {

        depth_val = gc->phong->Idist[(int)dist];

    }

    return (depth_val);
}


Vector ObjectNormalOnTheFly(GraphicalContext *gc, Image *image, int po, AdjRel *A)
{
    Vector N, V, valid_neighbor[A->n];
    int i, p, q, n, qo;
    Voxel u, v, uo, vo;

    p  = gc->surf_render[po].voxel;
    u  = GetVoxelCoord(gc->label, p);

    uo = GetVoxelCoord(image, po);

    n = 0;
    for (i = 1; i < A->n; i++)
    {
        vo = GetAdjacentVoxel(A, uo, i);
        if (ValidVoxel(image, vo))
        {
            qo = GetVoxelIndex(image, vo);
            q  = gc->surf_render[qo].voxel;
            if (q != NIL)
            {
                v  = GetVoxelCoord(gc->label, q);
                if (gc->label->val[q] == gc->label->val[p])
                {
                    valid_neighbor[n].x = v.x - u.x;
                    valid_neighbor[n].y = v.y - u.y;
                    valid_neighbor[n].z = v.z - u.z;
                    n++;
                }
            }
        }
    }

    N.x = N.y = N.z = 0.0;
    if ( n > 2 )
    {
        for (i = 0; i < n - 1; i++)
        {
            V = VectorCrossProduct(valid_neighbor[i], valid_neighbor[i + 1]);
            N.x += V.x;
            N.y += V.y;
            N.z += V.z;
        }
        N.x /= -n; N.y /= -n; N.z /= -n;
        N = NormalizeVector(N);
    }

    return (N);
}

void RenderFromSRBuffers(GraphicalContext *gc, Image *image)
{
    AdjRel *A = ClockCircular(3.0);

    for (int po = 0; po < image->n; po++)
    {

        Vector N;
        float     opac, dist, phong_val;
        int       p;
        Color  RGB, YCbCr;

        if (gc->surf_render[po].voxel != NIL)
        {

            N          = ObjectNormalOnTheFly(gc, image, po, A);
            opac       = gc->object[gc->surf_render[po].object].opacity;
            p          = gc->surf_render[po].voxel;
            dist       = gc->surf_render[po].depth;

            phong_val  = opac * PhongShading(gc, p, N, dist);
            RGB.val[0] = (int)(255.0 * phong_val * gc->object[gc->surf_render[po].object].red);
            RGB.val[1] = (int)(255.0 * phong_val * gc->object[gc->surf_render[po].object].green);
            RGB.val[2] = (int)(255.0 * phong_val * gc->object[gc->surf_render[po].object].blue);
            YCbCr      = RGBtoYCbCr(RGB);

            image->val[po] = YCbCr.val[0];
            image->Cb[po]  = YCbCr.val[1];
            image->Cr[po]  = YCbCr.val[2];

        }
    }

    DestroyAdjRel(&A);
}

void ResetSRBuffers(GraphicalContext *gc)
{
    int n = DiagonalSize(gc->label);

    n = n * n;
    for (int p = 0; p < n; p++)
    {
        gc->surf_render[p].depth    = INFINITY_FLT;
        gc->surf_render[p].opacity  = 1.0;
        gc->surf_render[p].voxel    = NIL;
        gc->surf_render[p].object   = NIL;
    }
}

float OpacityValueAtPoint(GraphicalContext *gc, Point P)
{
    Voxel u[8];
    int      p[8], i;
    float    dx = 1.0, dy = 1.0, dz = 1.0, val[6];
    float    value;
    float    o[8];

    if ((int)(P.x + 1.0) == gc->scene->xsize) dx = 0.0;
    if ((int)(P.y + 1.0) == gc->scene->ysize) dy = 0.0;
    if ((int)(P.z + 1.0) == gc->scene->zsize) dz = 0.0;

    u[0].x = (int)P.x;      u[0].y = (int)P.y;       u[0].z = (int)P.z;
    u[1].x = (int)(P.x + dx); u[1].y = (int)P.y;       u[1].z = (int)P.z;
    u[2].x = (int)P.x;      u[2].y = (int)(P.y + dy);  u[2].z = (int)P.z;
    u[3].x = (int)(P.x + dx); u[3].y = (int)(P.y + dy);  u[3].z = (int)P.z;
    u[4].x = (int)P.x;      u[4].y = (int)P.y;       u[4].z = (int)(P.z + dz);
    u[5].x = (int)(P.x + dx); u[5].y = (int)P.y;       u[5].z = (int)(P.z + dz);
    u[6].x = (int)P.x;      u[6].y = (int)(P.y + dy);  u[6].z = (int)(P.z + dz);
    u[7].x = (int)(P.x + dx); u[7].y = (int)(P.y + dy);  u[7].z = (int)(P.z + dz);

    if (gc->opacity != NULL)  /* volume rendering */
    {
        for (i = 0; i < 8; i++)
        {
            if (FValidVoxel(gc->scene, u[i]))
            {
                p[i] = FGetVoxelIndex(gc->scene, u[i]);
                o[i] = gc->opacity->val[p[i]];
            }
            else
            {
                u[0].x = ROUND(P.x);
                u[0].y = ROUND(P.y);
                u[0].z = ROUND(P.z);
                p[0]   = FGetVoxelIndex(gc->scene, u[0]);
                return (gc->opacity->val[p[0]]);
            }
        }
    }
    else   /* surface rendering */
    {
        if (gc->nobjs > 0)
        {
            for (i = 0; i < 8; i++)
            {
                if (FValidVoxel(gc->scene, u[i]))
                {
                    p[i] = FGetVoxelIndex(gc->scene, u[i]);
                    o[i] = gc->object[gc->label->val[p[i]]].opacity;
                }
                else
                {
                    u[0].x = ROUND(P.x);
                    u[0].y = ROUND(P.y);
                    u[0].z = ROUND(P.z);
                    p[0]   = FGetVoxelIndex(gc->scene, u[0]);
                    return (gc->object[gc->label->val[p[0]]].opacity);
                }
            }
        }
        else
        {
            Error("Cannot interpolate opacities", "OpacityValueAtPoint");
        }
    }

    val[0] = (float)o[1] * (P.x - u[0].x) + (float)o[0] * (u[1].x - P.x);
    val[1] = (float)o[3] * (P.x - u[2].x) + (float)o[2] * (u[3].x - P.x);
    val[2] = (float)o[5] * (P.x - u[4].x) + (float)o[4] * (u[5].x - P.x);
    val[3] = (float)o[7] * (P.x - u[6].x) + (float)o[6] * (u[7].x - P.x);
    val[4] = val[1] * (P.y - u[0].y) + val[0] * (u[2].y - P.y);
    val[5] = val[3] * (P.y - u[0].y) + val[2] * (u[2].y - P.y);
    value  = (val[5] * (P.z - u[0].z) + val[4] * (u[4].z - P.z));

    return (value);
}

Vector NormalVectorAtPoint(GraphicalContext *gc, Point P)
{
    Voxel  u[8];
    int       p[8], i;
    float     dx = 1.0, dy = 1.0, dz = 1.0, val[8];
    Vector V[8];

    if ((int)(P.x + 1.0) == gc->scene->xsize) dx = 0.0;
    if ((int)(P.y + 1.0) == gc->scene->ysize) dy = 0.0;
    if ((int)(P.z + 1.0) == gc->scene->zsize) dz = 0.0;

    u[0].x = (int)P.x;      u[0].y = (int)P.y;       u[0].z = (int)P.z;
    u[1].x = (int)(P.x + dx); u[1].y = (int)P.y;       u[1].z = (int)P.z;
    u[2].x = (int)P.x;      u[2].y = (int)(P.y + dy);  u[2].z = (int)P.z;
    u[3].x = (int)(P.x + dx); u[3].y = (int)(P.y + dy);  u[3].z = (int)P.z;
    u[4].x = (int)P.x;      u[4].y = (int)P.y;       u[4].z = (int)(P.z + dz);
    u[5].x = (int)(P.x + dx); u[5].y = (int)P.y;       u[5].z = (int)(P.z + dz);
    u[6].x = (int)P.x;      u[6].y = (int)(P.y + dy);  u[6].z = (int)(P.z + dz);
    u[7].x = (int)(P.x + dx); u[7].y = (int)(P.y + dy);  u[7].z = (int)(P.z + dz);

    for (i = 0; i < 8; i++)
    {
        if (FValidVoxel(gc->scene, u[i]))
        {
            p[i] = FGetVoxelIndex(gc->scene, u[i]);
            V[i].x = gc->phong->normal[gc->normal->val[p[i]]].x;
            V[i].y = gc->phong->normal[gc->normal->val[p[i]]].y;
            V[i].z = gc->phong->normal[gc->normal->val[p[i]]].z;
        }
        else
        {
            u[0].x = ROUND(P.x);
            u[0].y = ROUND(P.y);
            u[0].z = ROUND(P.z);
            p[0]   = FGetVoxelIndex(gc->scene, u[0]);
            V[0].x = gc->phong->normal[gc->normal->val[p[i]]].x;
            V[0].y = gc->phong->normal[gc->normal->val[p[i]]].y;
            V[0].z = gc->phong->normal[gc->normal->val[p[i]]].z;
            return (V[0]);
        }
    }

    val[0] = (float)V[1].x * (P.x - u[0].x) + (float)V[0].x * (u[1].x - P.x);
    val[1] = (float)V[3].x * (P.x - u[2].x) + (float)V[2].x * (u[3].x - P.x);
    val[2] = (float)V[5].x * (P.x - u[4].x) + (float)V[4].x * (u[5].x - P.x);
    val[3] = (float)V[7].x * (P.x - u[6].x) + (float)V[6].x * (u[7].x - P.x);
    val[4] = val[1] * (P.y - u[0].y) + val[0] * (u[2].y - P.y);
    val[5] = val[3] * (P.y - u[0].y) + val[2] * (u[2].y - P.y);
    V[0].x = (val[5] * (P.z - u[0].z) + val[4] * (u[4].z - P.z));

    val[0] = (float)V[1].y * (P.x - u[0].x) + (float)V[0].y * (u[1].x - P.x);
    val[1] = (float)V[3].y * (P.x - u[2].x) + (float)V[2].y * (u[3].x - P.x);
    val[2] = (float)V[5].y * (P.x - u[4].x) + (float)V[4].y * (u[5].x - P.x);
    val[3] = (float)V[7].y * (P.x - u[6].x) + (float)V[6].y * (u[7].x - P.x);
    val[4] = val[1] * (P.y - u[0].y) + val[0] * (u[2].y - P.y);
    val[5] = val[3] * (P.y - u[0].y) + val[2] * (u[2].y - P.y);
    V[0].y = (val[5] * (P.z - u[0].z) + val[4] * (u[4].z - P.z));

    val[0] = (float)V[1].z * (P.x - u[0].x) + (float)V[0].z * (u[1].x - P.x);
    val[1] = (float)V[3].z * (P.x - u[2].x) + (float)V[2].z * (u[3].x - P.x);
    val[2] = (float)V[5].z * (P.x - u[4].x) + (float)V[4].z * (u[5].x - P.x);
    val[3] = (float)V[7].z * (P.x - u[6].x) + (float)V[6].z * (u[7].x - P.x);
    val[4] = val[1] * (P.y - u[0].y) + val[0] * (u[2].y - P.y);
    val[5] = val[3] * (P.y - u[0].y) + val[2] * (u[2].y - P.y);
    V[0].z = (val[5] * (P.z - u[0].z) + val[4] * (u[4].z - P.z));

    return (V[0]);
}

void SurfaceShadingAlongRay(GraphicalContext *gc, Point P0, Point P1, Point Pn, float *red, float *green, float *blue)
{
    Point     u;
    int          k, n, p;
    float        Dx = (Pn.x - P1.x), Dy = (Pn.y - P1.y), Dz = (Pn.z - P1.z), dist;
    float        dx = 0, dy = 0, dz = 0, phong_val, acum_opacity, opac;
    Voxel     v;
    Vector    N;

    *red = *green = *blue = 0.0;

    /* DDA - Digital Differential Analyzer */

    if (PointsAreEqual(P1, Pn))
    {
        n = 1;
    }
    else   /* process points from P1 to Pn */
    {
        if ((fabs(Dx) >= fabs(Dy)) && (fabs(Dx) >= fabs(Dz)))
        {
            /* Dx is the maximum projection of
                                         vector P1Pn */
            n  = (int)(fabs(Dx) + 1);
            dx = SIGN(Dx);
            dy = dx * Dy / Dx;
            dz = dx * Dz / Dx;
        }
        else
        {
            if ((fabs(Dy) >= fabs(Dx)) && (fabs(Dy) >= fabs(Dz)))
            {
                /* Dy is the maximum projection of
                                         vector P1Pn */
                n  = (int)(fabs(Dy) + 1);
                dy = SIGN(Dy);
                dx = dy * Dx / Dy;
                dz = dy * Dz / Dy;
            }
            else     /* Dz is the maximum projection of vector P1Pn */
            {
                n  = (int)(fabs(Dz) + 1);
                dz = SIGN(Dz);
                dx = dz * Dx / Dz;
                dy = dz * Dy / Dz;
            }
        }
    }

    /* Execute shading model along the viewing ray */

    u.x  = P1.x;  u.y = P1.y;  u.z = P1.z;

    for (k = 0, acum_opacity = 1.0; (k < n) && (acum_opacity > Epsilon); k++)
    {

        v.x = ROUND(u.x);   v.y = ROUND(u.y);   v.z = ROUND(u.z);
        p = FGetVoxelIndex(gc->scene, v);
        if (gc->label->val[p] != 0)
        {
            if (gc->object[gc->label->val[p]].visibility != 0)
            {
                if (gc->object[gc->label->val[p]].opacity > 0.0)
                {
                    dist = PointDistance(u, P0);
                    N.x  = gc->phong->normal[gc->normal->val[p]].x;
                    N.y  = gc->phong->normal[gc->normal->val[p]].y;
                    N.z  = gc->phong->normal[gc->normal->val[p]].z;
                    opac = gc->object[gc->label->val[p]].opacity;
                    phong_val = opac * PhongShading(gc, p, N, dist) * acum_opacity;
                    *red   += phong_val * gc->object[gc->label->val[p]].red;
                    *green += phong_val * gc->object[gc->label->val[p]].green;
                    *blue  += phong_val * gc->object[gc->label->val[p]].blue;
                    acum_opacity = acum_opacity * (1.0 -  opac);
                }
            }
        }
        u.x = u.x + dx;
        u.y = u.y + dy;
        u.z = u.z + dz;
    }


    if (*red < 0) *red = 0;
    if (*green < 0) *green = 0;
    if (*blue < 0) *blue = 0;
    if (*red > 1) *red = 1;
    if (*green > 1) *green = 1;
    if (*blue > 1) *blue = 1;

}

void VoxelInfoByRayCasting(GraphicalContext *gc, Point P0, Point P1, Point Pn, int *voxel, float *depth, int *object) /* only for scenes with opaque objects */
{
    Point     u;
    int          k, n, p;
    float        Dx = (Pn.x - P1.x), Dy = (Pn.y - P1.y), Dz = (Pn.z - P1.z);
    float        dx = 0, dy = 0, dz = 0;
    Voxel     v;

    *voxel = NIL; *depth = INFINITY_FLT; *object = NIL;

    /* DDA - Digital Differential Analyzer */

    if (PointsAreEqual(P1, Pn))
    {
        n = 1;
    }
    else   /* process points from P1 to Pn */
    {

        if ((fabs(Dx) >= fabs(Dy)) && (fabs(Dx) >= fabs(Dz)))
        {
            /* Dx is the maximum projection of
                                         vector P1Pn */
            n  = (int)(fabs(Dx) + 1);
            dx = SIGN(Dx);
            dy = dx * Dy / Dx;
            dz = dx * Dz / Dx;
        }
        else
        {
            if ((fabs(Dy) >= fabs(Dx)) && (fabs(Dy) >= fabs(Dz)))
            {
                /* Dy is the maximum projection of
                                         vector P1Pn */
                n  = (int)(fabs(Dy) + 1);
                dy = SIGN(Dy);
                dx = dy * Dx / Dy;
                dz = dy * Dz / Dy;
            }
            else     /* Dz is the maximum projection of vector P1Pn */
            {
                n  = (int)(fabs(Dz) + 1);
                dz = SIGN(Dz);
                dx = dz * Dx / Dz;
                dy = dz * Dy / Dz;
            }
        }
    }


    /* Find the closest voxel and its depth to the viewing plane */

    u.x  = P1.x;  u.y = P1.y;  u.z = P1.z;


    for (k = 0; k < n; k++)
    {
        v.x = ROUND(u.x);   v.y = ROUND(u.y);   v.z = ROUND(u.z);
        p = FGetVoxelIndex(gc->scene, v);
        if (gc->label->val[p] != 0)
        {
            if (gc->object[gc->label->val[p]].visibility != 0)
            {
                if (gc->object[gc->label->val[p]].opacity > 0.0)
                {
                    *depth  = PointDistance(u, P0);
                    *voxel  = p;
                    *object = gc->label->val[p];
                    break;
                }
            }
        }

        u.x = u.x + dx;
        u.y = u.y + dy;
        u.z = u.z + dz;
    }

}


void VolumeShadingAlongRay(GraphicalContext *gc, Point P0, Point P1, Point Pn, float *value)
{
    Point     u;
    int          k, n, p;
    float        Dx = (Pn.x - P1.x), Dy = (Pn.y - P1.y), Dz = (Pn.z - P1.z);
    float        dx = 0, dy = 0, dz = 0, acum_opacity = 1.0, dist, opac;
    Vector    N;
    Voxel     v;

    *value = 0.0;

    /* DDA - Digital Differential Analyzer */

    if (PointsAreEqual(P1, Pn))
    {
        n = 1;
    }
    else   /* process points from P1 to Pn */
    {
        if ((fabs(Dx) >= fabs(Dy)) && (fabs(Dx) >= fabs(Dz)))
        {
            /* Dx is the maximum projection of
                                         vector P1Pn */
            n  = (int)(fabs(Dx) + 1);
            dx = SIGN(Dx);
            dy = dx * Dy / Dx;
            dz = dx * Dz / Dx;
        }
        else
        {
            if ((fabs(Dy) >= fabs(Dx)) && (fabs(Dy) >= fabs(Dz)))
            {
                /* Dy is the maximum projection of
                                         vector P1Pn */
                n  = (int)(fabs(Dy) + 1);
                dy = SIGN(Dy);
                dx = dy * Dx / Dy;
                dz = dy * Dz / Dy;
            }
            else     /* Dz is the maximum projection of vector P1Pn */
            {
                n  = (int)(fabs(Dz) + 1);
                dz = SIGN(Dz);
                dx = dz * Dx / Dz;
                dy = dz * Dy / Dz;
            }
        }
    }

    u.x  = P1.x;  u.y = P1.y;  u.z = P1.z;

    /* Execute Phong model along the viewing ray */

    for (k = 0; (k < n) && (acum_opacity > Epsilon); k++)
    {
        v.x = ROUND(u.x);   v.y = ROUND(u.y);   v.z = ROUND(u.z);
        p = FGetVoxelIndex(gc->scene, v);
        if (gc->opacity->val[p] > 0.0)
        {
            dist = PointDistance(u, P0);
            N.x  = gc->phong->normal[gc->normal->val[p]].x;
            N.y  = gc->phong->normal[gc->normal->val[p]].y;
            N.z  = gc->phong->normal[gc->normal->val[p]].z;
            N    = NormalVectorAtPoint(gc, u);
            opac = OpacityValueAtPoint(gc, u);
            *value += opac * PhongShading(gc, p, N, dist) * acum_opacity;
            acum_opacity = acum_opacity * (1.0 -  opac);
        }
        u.x = u.x + dx;
        u.y = u.y + dy;
        u.z = u.z + dz;
    }

}


Image *SurfaceRenderingByRayCasting(GraphicalContext *gc)
{
    Image  *image;
    Vector  n;
    int        diag  = FDiagonalSize(gc->scene);

    n.x  = 0; n.y = 0; n.z = 1;
    n    = TransformVector(gc->viewdir->Rinv, n);

    image =  CreateImage(diag, diag, 1);
    SetCbCr(image, 128);

    if (gc->normal != NULL)
    {

        #pragma omp parallel for shared(image,gc,n)
        for (int po = 0; po < image->n; po++)
        {

            Point  P0, P1, Pn;
            Color  RGB, YCbCr;
            float     red = 0.0, green = 0.0, blue = 0.0;

            P0.x   =  GetXCoord(image, po);
            P0.y   =  GetYCoord(image, po);
            P0.z   = -diag / 2.0;
            P0     =  TransformPoint(gc->viewdir->Tinv, P0);

            if (IntersectionPoints(gc, P0, n, &P1, &Pn))
            {

                SurfaceShadingAlongRay(gc, P0, P1, Pn, &red, &green, &blue);

                RGB.val[0]     = (int)(255.0 * red);
                RGB.val[1]     = (int)(255.0 * green);
                RGB.val[2]     = (int)(255.0 * blue);
                YCbCr          = RGBtoYCbCr(RGB);
                image->val[po] = YCbCr.val[0];
                image->Cb[po]  = YCbCr.val[1];
                image->Cr[po]  = YCbCr.val[2];

            }
        }
    }
    else     /* Estimate normals on-the-fly */
    {

        if (gc->overall_opac == 1.0)
        {

            #pragma omp parallel for shared(image,gc,n)
            for (int po = 0; po < image->n; po++)
            {

                Point  P0, P1, Pn;

                P0.x   =  GetXCoord(image, po);
                P0.y   =  GetYCoord(image, po);
                P0.z   = -diag / 2.0;
                P0     =  TransformPoint(gc->viewdir->Tinv, P0);

                if (IntersectionPoints(gc, P0, n, &P1, &Pn))
                {

                    VoxelInfoByRayCasting(gc, P0, P1, Pn, &gc->surf_render[po].voxel, &gc->surf_render[po].depth, &gc->surf_render[po].object);


                }
            }

            RenderFromSRBuffers(gc, image);

        }
        else
        {
            Error("Semi-transparent objects require to set object normals", "SurfaceRenderingByRayCasting");
        }

    }

    return (image);
}

Image *SurfaceRenderingBySplatting(GraphicalContext *gc)
{
    int        xo, yo, zo, xf, yf, zf, dx, dy, dz, p, diag = FDiagonalSize(gc->scene);
    Voxel   u;
    Point   P;
    Vector  N;
    float      opac, phong_val;
    Image  *image = CreateImage(diag, diag, 1);

    SetCbCr(image, 128);

    xo = gc->viewdir->ftb[gc->viewdir->octant].xo; yo = gc->viewdir->ftb[gc->viewdir->octant].yo; zo = gc->viewdir->ftb[gc->viewdir->octant].zo;
    xf = gc->viewdir->ftb[gc->viewdir->octant].xf; yf = gc->viewdir->ftb[gc->viewdir->octant].yf; zf = gc->viewdir->ftb[gc->viewdir->octant].zf;
    dx = gc->viewdir->ftb[gc->viewdir->octant].dx; dy = gc->viewdir->ftb[gc->viewdir->octant].dy; dz = gc->viewdir->ftb[gc->viewdir->octant].dz;

    switch (gc->viewdir->paxis)
    {

    case AXIS_X:

        for (u.x = xo; (u.x != xf); u.x = u.x + dx)
            for (u.y = yo; (u.y != yf); u.y = u.y + dy)
                for (u.z = zo; (u.z != zf); u.z = u.z + dz)
                {
                    p = FGetVoxelIndex(gc->scene, u);
                    if ((gc->opacity->val[p] > 0.0) &&
                            (gc->object[gc->label->val[p]].opacity > 0.0) &&
                            (gc->object[gc->label->val[p]].visibility != 0))
                    {
                        P.x       = u.x; P.y = u.y; P.z = u.z;
                        P         = TransformPoint(gc->viewdir->T, P);
                        N.x       = gc->phong->normal[gc->normal->val[p]].x;
                        N.y       = gc->phong->normal[gc->normal->val[p]].y;
                        N.z       = gc->phong->normal[gc->normal->val[p]].z;
                        N         = TransformVector(gc->viewdir->R, N);
                        opac      = gc->object[gc->label->val[p]].opacity;
                        phong_val = opac * PhongShading(gc, p, N, P.z);
                        VoxelSplatting(gc, image, P, p, phong_val, opac);
                    }
                }

        break;

    case AXIS_Y:

        for (u.y = yo; (u.y != yf); u.y = u.y + dy)
            for (u.x = xo; (u.x != xf); u.x = u.x + dx)
                for (u.z = zo; (u.z != zf); u.z = u.z + dz)
                {
                    p = FGetVoxelIndex(gc->scene, u);
                    if ((gc->opacity->val[p] != 0) &&
                            (gc->object[gc->label->val[p]].opacity > 0.0) &&
                            (gc->object[gc->label->val[p]].visibility != 0))
                    {
                        P.x       = u.x; P.y = u.y; P.z = u.z;
                        P         = TransformPoint(gc->viewdir->T, P);
                        N.x       = gc->phong->normal[gc->normal->val[p]].x;
                        N.y       = gc->phong->normal[gc->normal->val[p]].y;
                        N.z       = gc->phong->normal[gc->normal->val[p]].z;
                        N         = TransformVector(gc->viewdir->R, N);
                        opac      = gc->object[gc->label->val[p]].opacity;
                        phong_val = opac * PhongShading(gc, p, N, P.z);
                        VoxelSplatting(gc, image, P, p, phong_val, opac);
                    }
                }

        break;

    case AXIS_Z:

        for (u.z = zo; (u.z != zf); u.z = u.z + dz)
            for (u.x = xo; (u.x != xf); u.x = u.x + dx)
                for (u.y = yo; (u.y != yf); u.y = u.y + dy)
                {
                    p = FGetVoxelIndex(gc->scene, u);
                    if ((gc->opacity->val[p] != 0) &&
                            (gc->object[gc->label->val[p]].opacity > 0.0) &&
                            (gc->object[gc->label->val[p]].visibility != 0))
                    {
                        P.x       = u.x; P.y = u.y; P.z = u.z;
                        P         = TransformPoint(gc->viewdir->T, P);
                        N.x       = gc->phong->normal[gc->normal->val[p]].x;
                        N.y       = gc->phong->normal[gc->normal->val[p]].y;
                        N.z       = gc->phong->normal[gc->normal->val[p]].z;
                        N         = TransformVector(gc->viewdir->R, N);
                        opac      = gc->object[gc->label->val[p]].opacity;
                        phong_val = opac * PhongShading(gc, p, N, P.z);
                        VoxelSplatting(gc, image, P, p, phong_val, opac);
                    }
                }
        break;

    default:
        Error("Invalid principal axis", "SurfaceRenderingBySplatting");

    }

    return (image);
}

Image *VolumeRenderingByRayCasting(GraphicalContext *gc)
{
    Image  *image;
    Vector  n;
    int        diag  = FDiagonalSize(gc->scene);

    image =  CreateImage(diag, diag, 1);

    n.x  = 0; n.y = 0; n.z = 1;
    n    = TransformVector(gc->viewdir->Rinv, n);

    for (int po = 0; po < image->n; po++)
    {
        Point  P0, P1, Pn;
        float     value;

        P0.x   =  GetXCoord(image, po);
        P0.y   =  GetYCoord(image, po);
        P0.z   = -diag / 2.0;
        P0     =  TransformPoint(gc->viewdir->Tinv, P0);

        if (IntersectionPoints(gc, P0, n, &P1, &Pn))
        {
            VolumeShadingAlongRay(gc, P0, P1, Pn, &value);
            image->val[po] = (int)(255.0 * value);
        }
    }

    return (image);
}

void SetSceneFaces(GraphicalContext *gc)
{
    FImage *scene = gc->scene;

    gc->face[0].pos.x = scene->xsize / 2.0;
    gc->face[0].pos.y = scene->ysize / 2.0;
    gc->face[0].pos.z = 0;
    gc->face[1].pos.x = scene->xsize / 2.0;
    gc->face[1].pos.y = scene->ysize / 2.0;
    gc->face[1].pos.z = scene->zsize - 1;
    gc->face[2].pos.x = scene->xsize / 2.0;
    gc->face[2].pos.y = 0;
    gc->face[2].pos.z = scene->zsize / 2.0;
    gc->face[3].pos.x = scene->xsize / 2.0;
    gc->face[3].pos.y = scene->ysize - 1;
    gc->face[3].pos.z = scene->zsize / 2.0;
    gc->face[4].pos.x = 0;
    gc->face[4].pos.y = scene->ysize / 2.0;
    gc->face[4].pos.z = scene->zsize / 2.0;
    gc->face[5].pos.x = scene->xsize - 1;
    gc->face[5].pos.y = scene->ysize / 2.0;
    gc->face[5].pos.z = scene->zsize / 2.0;

    gc->face[0].normal.x =  0;
    gc->face[0].normal.y =  0;
    gc->face[0].normal.z = -1;
    gc->face[1].normal.x =  0;
    gc->face[1].normal.y =  0;
    gc->face[1].normal.z =  1;
    gc->face[2].normal.x =  0;
    gc->face[2].normal.y = -1;
    gc->face[2].normal.z =  0;
    gc->face[3].normal.x =  0;
    gc->face[3].normal.y =  1;
    gc->face[3].normal.z =  0;
    gc->face[4].normal.x = -1;
    gc->face[4].normal.y =  0;
    gc->face[4].normal.z =  0;
    gc->face[5].normal.x =  1;
    gc->face[5].normal.y =  0;
    gc->face[5].normal.z =  0;

}

Vector *CreateNormalTable()
{
    int        i, gamma, alpha;
    float      gamma_rad, alpha_rad;
    Vector *normaltable;

    /* creates normal look-up table */

    normaltable = (Vector *)calloc(NUM_OF_NORMALS, sizeof(Vector));

    if (normaltable == NULL)
    {
        Error(MSG1, "CreateNormalTable");
    }

    normaltable[0].x = 0.0;
    normaltable[0].y = 0.0;
    normaltable[0].z = 0.0;

    i = 1;
    for (gamma = -90; gamma <= 90; gamma++)
    {
        gamma_rad = (PI * gamma) / 180.0;
        for (alpha = 0; alpha < 360; alpha++)
        {
            alpha_rad = (PI * alpha) / 180.0;
            normaltable[i].x = cos(gamma_rad) * cos(alpha_rad);
            normaltable[i].y = cos(gamma_rad) * sin(alpha_rad);
            normaltable[i].z = sin(gamma_rad);
            i++;
        }
    }

    return (normaltable);
}

int GetNormalIndex(Vector normal)
{
    int gamma, alpha, index;

    if ((normal.x == 0.0) && (normal.y == 0.0) && (normal.z == 0.0))
    {
        return (0);
    }

    gamma = (int)(asin(normal.z) * 180.0 / PI); /* [-90,90] */
    alpha = (int)(atan2(normal.y, normal.x) * 180.0 / PI); /* [-180,180] */
    if (alpha < 0)
        alpha += 360;
    index = ((gamma + 90) * 360) + alpha + 1;

    return (index);
}

ObjectAttributes *CreateObjectAttributes(Image *label, int *nobjs)
{
    ObjectAttributes *object;
    *nobjs = MaximumValue(label);

    object = (ObjectAttributes *)calloc(*nobjs + 1, sizeof(ObjectAttributes));


    /* background */

    object[0].opacity    = 0;
    object[0].red        = 0;
    object[0].green      = 0;
    object[0].blue       = 0;
    object[0].visibility = 0;

    /* default for objects */

    for (int i = 1; i <= *nobjs; i++)
    {
        object[i].opacity    = 1;
        object[i].red        = 1;
        object[i].green      = 1;
        object[i].blue       = 1;
        object[i].visibility = 1;
    }

    return (object);
}

SRBuffers *CreateSRBuffers(Image *label)
{
    SRBuffers *surf_render;
    int n = DiagonalSize(label);

    n = n * n;
    surf_render            = (SRBuffers *) calloc(n, sizeof(SRBuffers));

    for (int p = 0; p < n; p++)
    {
        surf_render[p].depth    = INFINITY_FLT;
        surf_render[p].opacity  = 1.0;
        surf_render[p].voxel    = NIL;
        surf_render[p].object   = NIL;
    }

    return (surf_render);
}

PhongModel *CreatePhongModel(FImage *scene)
{
    PhongModel *phong = (PhongModel *)calloc(1, sizeof(PhongModel));

    phong->ka     = 0.1;
    phong->kd     = 0.7;
    phong->ks     = 0.2;
    phong->ns     = 5.0;
    phong->normal = CreateNormalTable();
    phong->ndists = (int)(2.0 * FDiagonalSize(scene) + 1);
    phong->Idist  = AllocFloatArray(phong->ndists);

    for (int d = 0; d < phong->ndists; d++)
    {
        phong->Idist[d] = (float)0.8 * (phong->ndists - d) / (float)phong->ndists + 0.2;
    }

    return (phong);
}

void SetFTBaccess(GraphicalContext *gc)
{
    gc->viewdir->ftb    = (FTBaccess *)calloc(8, sizeof(FTBaccess));

    /* Set the front-to-back voxel access from the closest to the farthest octant */

    int dx[8], dy[8], dz[8];

    dx[0] = 1; dx[1] = -1; dx[2] =  1; dx[3] = -1; dx[4] =  1; dx[5] = -1; dx[6] =  1; dx[7] = -1;
    dy[0] = 1; dy[1] =  1; dy[2] = -1; dy[3] = -1; dy[4] =  1; dy[5] =  1; dy[6] = -1; dy[7] = -1;
    dz[0] = 1; dz[1] =  1; dz[2] =  1; dz[3] =  1; dz[4] = -1; dz[5] = -1; dz[6] = -1; dz[7] = -1;

    int xo[8], yo[8], zo[8];

    xo[0] = 0; xo[1] = gc->scene->xsize - 1; xo[2] = 0;                  xo[3] = gc->scene->xsize - 1; xo[4] = 0;                  xo[5] = gc->scene->xsize - 1; xo[6] = 0;                  xo[7] = gc->scene->xsize - 1;
    yo[0] = 0; yo[1] = 0;                  yo[2] = gc->scene->ysize - 1; yo[3] = gc->scene->ysize - 1; yo[4] = 0;                  yo[5] = 0;                  yo[6] = gc->scene->ysize - 1; yo[7] = gc->scene->ysize - 1;
    zo[0] = 0; zo[1] = 0;                  zo[2] = 0;                  zo[3] = 0;                  zo[4] = gc->scene->zsize - 1; zo[5] = gc->scene->zsize - 1; zo[6] = gc->scene->zsize - 1; zo[7] = gc->scene->zsize - 1;

    int xf[8], yf[8], zf[8];

    for (int i = 0; i < 8; i++)
    {
        if (dx[i] == 1)
            xf[i] = gc->scene->xsize;
        else
            xf[i] = -1;
        if (dy[i] == 1)
            yf[i] = gc->scene->ysize;
        else
            yf[i] = -1;
        if (dz[i] == 1)
            zf[i] = gc->scene->zsize;
        else
            zf[i] = -1;
    }

    for (int i = 0; i < 8; i++)
    {
        gc->viewdir->ftb[i].dx = dx[i];    gc->viewdir->ftb[i].dy = dy[i];    gc->viewdir->ftb[i].dz = dz[i];
        gc->viewdir->ftb[i].xo = xo[i];    gc->viewdir->ftb[i].yo = yo[i];    gc->viewdir->ftb[i].zo = zo[i];
        gc->viewdir->ftb[i].xf = xf[i];    gc->viewdir->ftb[i].yf = yf[i];    gc->viewdir->ftb[i].zf = zf[i];
    }

}

/* ------------------------------- Public Methods --------------------------*/

void SetSceneNormal(GraphicalContext *gc)
{
    AdjRel *A   = Spheric(3.0);
    float     *mag = AllocFloatArray(A->n);
    float      Delta;
    int        i, p, q;
    Voxel   u, v;
    Vector  N;


    if (gc->opacity == NULL)
        Error("Set scene opacity first", "SetSceneNormal");

    if (gc->normal != NULL)
        DestroyImage(&gc->normal);

    gc->normal = CreateImage(gc->scene->xsize, gc->scene->ysize, gc->scene->zsize);

    for (i = 0; i < A->n; i++)
        mag[i] = sqrtf(A->adj[i].dx * A->adj[i].dx + A->adj[i].dy * A->adj[i].dy + A->adj[i].dz * A->adj[i].dz);

    for (u.z = 0; u.z < gc->scene->zsize; u.z++)
        for (u.y = 0; u.y < gc->scene->ysize; u.y++)
            for (u.x = 0; u.x < gc->scene->xsize; u.x++)
            {
                p = FGetVoxelIndex(gc->scene, u);
                if (gc->opacity->val[p] > 0.0)
                {
                    N.x = N.y = N.z = 0.0;
                    for (i = 1; i < A->n; i++)
                    {
                        v = GetAdjacentVoxel(A, u, i);
                        if (FValidVoxel(gc->scene, v))
                        {
                            q = FGetVoxelIndex(gc->scene, v);
                            Delta = gc->scene->val[q] - gc->scene->val[p];
                            N.x  += Delta * A->adj[i].dx / mag[i];
                            N.y  += Delta * A->adj[i].dy / mag[i];
                            N.z  += Delta * A->adj[i].dz / mag[i];
                        }
                    }
                    /* it assumes scenes with objects brighter than the background */
                    N = NormalizeVector(N);
                    N.x = -N.x; N.y = -N.y; N.z = -N.z;
                    gc->normal->val[p] = GetNormalIndex(N);
                }
            }


    free(mag);
    DestroyAdjRel(&A);

}

void SetProjectionMode(GraphicalContext *gc, char proj_mode)
{
    gc->proj_mode = proj_mode;
    gc->viewdir->V.x = 0; gc->viewdir->V.y = 0; gc->viewdir->V.z = -1;
    if (gc->proj_mode == RAYCASTING)
    {
        gc->viewdir->V   = TransformVector(gc->viewdir->Rinv, gc->viewdir->V);
    }
}

void SetObjectNormal(GraphicalContext *gc)
{
    AdjRel *A;
    FImage *dist;
    float     *mag;
    float      Delta;
    int        i, p, q;
    Voxel   u, v;
    Vector  N;
    Set    *S = NULL;

    if (gc->label == NULL)
        Error("Object labels are required", "SetObjectNormal");

    if (gc->normal != NULL)
        DestroyImage(&gc->normal);

    /* extract shell around the border and some border voxels */

    A             = Spheric(sqrtf(3.0));
    dist          = ShellSignedDistTrans(gc->label, A, 5);
    S             = ObjectBorderSet(gc->label, A);
    DestroyAdjRel(&A);

    /* estimate object-based normal vectors and set opacity scene for the shell */

    gc->normal  = CreateImage(dist->xsize, dist->ysize, dist->zsize);
    gc->opacity = CreateFImage(dist->xsize, dist->ysize, dist->zsize);

    A   = Spheric(5.0);
    mag = AllocFloatArray(A->n);
    for (i = 0; i < A->n; i++)
        mag[i] = sqrtf(A->adj[i].dx * A->adj[i].dx + A->adj[i].dy * A->adj[i].dy + A->adj[i].dz * A->adj[i].dz);

    while (S != NULL)
    {

        p = RemoveSet(&S);

        gc->opacity->val[p] = 1.0;

        u = FGetVoxelCoord(dist, p);
        N.x = N.y = N.z = 0.0;

        for (i = 1; i < A->n; i++)
        {
            v = GetAdjacentVoxel(A, u, i);
            if (FValidVoxel(dist, v))
            {
                q = FGetVoxelIndex(dist, v);
                Delta = dist->val[q] - dist->val[p];
                N.x  += Delta * A->adj[i].dx / mag[i];
                N.y  += Delta * A->adj[i].dy / mag[i];
                N.z  += Delta * A->adj[i].dz / mag[i];
            }
        }

        /* force normal to point outward the object */
        N = NormalizeVector(N);
        v.x = ROUND(u.x + N.x); v.y = ROUND(u.y + N.y); v.z = ROUND(u.z + N.z);
        if (FValidVoxel(dist, v))
        {
            q = FGetVoxelIndex(dist, v);
            if (gc->label->val[q] != 0)
            {
                N.x = -N.x; N.y = -N.y; N.z = -N.z;
            }
        }
        gc->normal->val[p] = GetNormalIndex(N);
    }

    free(mag);
    DestroyAdjRel(&A);
    DestroyFImage(&dist);

}

void SetObjectColor(GraphicalContext *gc, int object, float red, float green, float blue)
{

    if ((object > gc->nobjs) || (object <= 0))
        Error("Invalid object", "SetObjectColor");

    if ((red < 0) || (red > 1) || (green < 0) || (green > 1) || (blue < 0) || (blue > 1))
        Error("Invalid color value", "SetObjectColor");

    gc->object[object].red   = red;
    gc->object[object].green = green;
    gc->object[object].blue  = blue;

}

void SetObjectOpacity(GraphicalContext *gc, int object, float opacity)
{

    if ((object > gc->nobjs) || (object <= 0))
        Error("Invalid object", "SetObjectOpacity");

    if ((opacity < 0) || (opacity > 1))
        Error("Invalid opacity value", "SetObjectOpacity");

    gc->object[object].opacity   = opacity;

    gc->overall_opac = 1.0;
    for (int i = 1; i <= gc->nobjs; i++)
        gc->overall_opac *=  gc->object[i].opacity;

}

void SetObjectVisibility(GraphicalContext *gc, int object, char visibility)
{

    if ((object > gc->nobjs) || (object <= 0))
        Error("Invalid object", "SetObjectVisibility");

    if ((visibility != 0) && (visibility != 1))
        Error("Invalid visibility value", "SetObjectVisibility");

    gc->object[object].visibility = visibility;

}

void SetViewDir(GraphicalContext *gc, float tilt, float spin)
{

    Vector  t;
    Matrix *Rx, *Ry, *Txyz, *Tuv, *aux;
    int        diag = FDiagonalSize(gc->scene);

    if (gc->viewdir == NULL)
    {
        gc->viewdir = (ViewDir *)calloc(1, sizeof(ViewDir));
        SetFTBaccess(gc);
    }
    else
    {
        DestroyMatrix(&gc->viewdir->T);
        DestroyMatrix(&gc->viewdir->Tinv);
        DestroyMatrix(&gc->viewdir->R);
        DestroyMatrix(&gc->viewdir->Rinv);
    }

    /* Set scene transformation and rotation matrices */

    t.x = -gc->scene->xsize / 2.0; t.y = -gc->scene->ysize / 2.0; t.z = -gc->scene->zsize / 2.0;
    Txyz            =  TranslationMatrix(t.x, t.y, t.z);
    Rx              =  RotationMatrix(AXIS_X, tilt);
    Ry              =  RotationMatrix(AXIS_Y, spin);
    gc->viewdir->R  =  MultMatrices(Ry, Rx);
    t.x =  diag / 2.0; t.y = diag / 2.0; t.z = diag / 2.0;
    Tuv             =  TranslationMatrix(t.x, t.y, t.z);

    aux             =  MultMatrices(gc->viewdir->R, Txyz);
    gc->viewdir->T  =  MultMatrices(Tuv, aux);

    DestroyMatrix(aux);
    DestroyMatrix(Txyz);
    DestroyMatrix(Tuv);
    DestroyMatrix(Rx);
    DestroyMatrix(Ry);

    /* Set scene transformation and rotation inverse matrices */

    t.x      = -diag / 2.0; t.y = -diag / 2.0; t.z = -diag / 2.0;
    Tuv      =  TranslationMatrix(t.x, t.y, t.z);
    t.x      =  gc->scene->xsize / 2.0; t.y = gc->scene->ysize / 2.0; t.z = gc->scene->zsize / 2.0;
    Txyz     =  TranslationMatrix(t.x, t.y, t.z);
    Rx       =  RotationMatrix(AXIS_X, -tilt);
    Ry       =  RotationMatrix(AXIS_Y, -spin);
    gc->viewdir->Rinv  =  MultMatrices(Rx, Ry);

    aux               =  MultMatrices(gc->viewdir->Rinv, Tuv);
    gc->viewdir->Tinv =  MultMatrices(Txyz, aux);

    DestroyMatrix(aux);
    DestroyMatrix(Txyz);
    DestroyMatrix(Tuv);
    DestroyMatrix(Rx);
    DestroyMatrix(Ry);

    /* Set viewing direction */

    gc->viewdir->V.x = 0; gc->viewdir->V.y = 0; gc->viewdir->V.z = -1;
    if (gc->proj_mode == RAYCASTING)
    {
        gc->viewdir->V   = TransformVector(gc->viewdir->Rinv, gc->viewdir->V);
    }

    /* Set principal axis */

    float max_proj = -INFINITY_FLT;
    Vector N;

    N.x = 0; N.y = 0; N.z = 1;
    N = TransformVector(gc->viewdir->Rinv, N);
    if (fabs(N.x) > max_proj)
    {
        gc->viewdir->paxis  = AXIS_X;
        max_proj            = fabs(N.x);
    }
    if (fabs(N.y) > max_proj)
    {
        gc->viewdir->paxis  = AXIS_Y;
        max_proj            = fabs(N.y);
    }
    if (fabs(N.z) > max_proj)
    {
        gc->viewdir->paxis  = AXIS_Z;
        max_proj            = fabs(N.z);
    }

    /* Set closest octant */

    float depth = INFINITY_FLT;

    for (int i = 0; i < 8; i++)
    {
        Point P;

        P.x = gc->viewdir->ftb[i].xo; P.y = gc->viewdir->ftb[i].yo; P.z = gc->viewdir->ftb[i].zo;
        P   = TransformPoint(gc->viewdir->T, P);
        if (P.z < depth)
        {
            gc->viewdir->octant = i;
            depth = P.z;
        }
    }

}

void SetSceneOpacity(GraphicalContext *gc, float min_val, float max_val, Image *grad, int grad_thres, float max_opac)
{
    if (gc->opacity == NULL)
        gc->opacity = CreateFImage(gc->scene->xsize, gc->scene->ysize, gc->scene->zsize);

    for (int p = 0; p < grad->n; p++)
    {
        if ((gc->scene->val[p] >= min_val) && (gc->scene->val[p] <= max_val))
            gc->opacity->val[p] = max_opac / (1.0 + exp(-grad->val[p] + grad_thres));
    }

}

char IntersectionPoints(GraphicalContext *gc, Point P0, Vector n, Point *P1, Point *Pn)
{
    Vector V;
    Point  P;
    int i;
    float lambda, a, b, lambda_1 = INFINITY_FLT, lambda_n = -INFINITY_FLT;

    P1->x = Pn->x = 0;
    P1->y = Pn->y = 0;
    P1->z = Pn->z = 0;

    for (i = 0; i < 6; i++)
    {
        b = VectorInnerProduct(gc->face[i].normal, n);
        if (fabs(b) > Epsilon)
        {
            V.x    = gc->face[i].pos.x - P0.x;
            V.y    = gc->face[i].pos.y - P0.y;
            V.z    = gc->face[i].pos.z - P0.z;
            a      = VectorInnerProduct(gc->face[i].normal, V);
            lambda = a / b;
            P.x = ROUND(P0.x + lambda * n.x);
            P.y = ROUND(P0.y + lambda * n.y);
            P.z = ROUND(P0.z + lambda * n.z);

            if (FValidPoint(gc->scene, P))
            {
                if (lambda < lambda_1)
                {
                    lambda_1 = lambda;
                }
                if (lambda > lambda_n)
                {
                    lambda_n = lambda;
                }
            }
        }
    }

    if (lambda_1 < lambda_n)
    {

        P1->x = ROUND(P0.x + lambda_1 * n.x);
        P1->y = ROUND(P0.y + lambda_1 * n.y);
        P1->z = ROUND(P0.z + lambda_1 * n.z);

        Pn->x = ROUND(P0.x + lambda_n * n.x);
        Pn->y = ROUND(P0.y + lambda_n * n.y);
        Pn->z = ROUND(P0.z + lambda_n * n.z);

        return (1);
    }
    else
    {
        return (0);
    }
}

GraphicalContext *CreateGraphicalContext(FImage *scene, Image *label)
{
    GraphicalContext *gc;

    gc = (GraphicalContext *) calloc(1, sizeof(GraphicalContext));

    gc->scene          = FCopyImage(scene);
    gc->phong          = CreatePhongModel(scene);
    gc->viewdir        = NULL;
    gc->proj_mode      = RAYCASTING;
    gc->nobjs          = 0;
    gc->overall_opac   = 1.0;
    gc->face           = (Plane *) malloc(sizeof(Plane) * 6);
    SetViewDir(gc, 0, 0);
    SetSceneFaces(gc);

    if (label != NULL)  /* for surface rendering */
    {
        gc->label       = CopyImage(label);
        gc->object      = CreateObjectAttributes(label, &gc->nobjs);
        gc->surf_render = CreateSRBuffers(label);
    }
    else   /* for volume rendering */
    {
        gc->label       = NULL;
        gc->object      = NULL;
        gc->surf_render = NULL;
    }

    gc->opacity       = NULL;
    gc->normal        = NULL;

    return (gc);
}

void DestroyGraphicalContext(GraphicalContext *gc)
{
    if (gc != NULL)
    {
        if (gc->nobjs != 0)
        {
            free(gc->object);
            free(gc->surf_render);
            DestroyImage(&gc->label);
        }
        if (gc->opacity != NULL)
            DestroyFImage(&gc->opacity);
        if (gc->normal != NULL)
            DestroyImage(&gc->normal);
        free(gc->phong->normal);
        free(gc->phong->Idist);
        free(gc->phong);
        free(gc->face);
        DestroyMatrix(&gc->viewdir->T);
        DestroyMatrix(&gc->viewdir->Tinv);
        DestroyMatrix(&gc->viewdir->R);
        DestroyMatrix(&gc->viewdir->Rinv);
        free(gc->viewdir->ftb);
        free(gc->viewdir);
        DestroyFImage(&gc->scene);
        free(gc);
    }
    else
    {
        Warning("Graphical context is already NULL", "DestroyGraphicalContext");
    }
}

Image *SurfaceRender(GraphicalContext *gc)
{
    Image  *image = NULL;

    if (gc->nobjs == 0)
        Error("There are no objects for visualization", "SurfaceRender");

    if (gc->normal == NULL)  /* normals estimated on the fly */
    {
        ResetSRBuffers(gc);
        if (gc->proj_mode == SPLATTING)
        {
            Warning("Changing projection mode to ray casting", "SurfaceRender");
            SetProjectionMode(gc, RAYCASTING);
        }
        image = SurfaceRenderingByRayCasting(gc);
    }
    else   /* pre-computed normals */
    {
        if (gc->proj_mode == SPLATTING)
            image = SurfaceRenderingBySplatting(gc);
        else
            image = SurfaceRenderingByRayCasting(gc);
    }

    return (image);
}

Image *VolumeRender(GraphicalContext *gc)
{
    Image  *image = NULL;

    if (gc->opacity == NULL)
    {
        Error("Set opacity scene", "VolumeRenderingByRayCasting");
    }

    if (gc->normal == NULL)
    {
        SetSceneNormal(gc);
    }

    if (gc->proj_mode == SPLATTING)
    {
        Warning("Changing projection mode to ray casting", "VolumeRender");
        SetProjectionMode(gc, RAYCASTING);
    }

    image = VolumeRenderingByRayCasting(gc);

    return (image);
}




// void DrawPoint(Image *img, Voxel u, Color YCbCr, AdjRel *B)
// {
//   int q,i;
//   Voxel v;

//   if (!ValidVoxel(img,u))
//     Error("Point is outside the image domain","DrawPoint");

//   if (img->Cb == NULL){
//     img->Cb = AllocUShortArray(img->n);
//     img->Cr = AllocUShortArray(img->n);
//     for (q=0; q < img->n; q++) {
//       img->Cb[q] = 128;
//       img->Cr[q] = 128;
//     }
//   }

//   for (i=0; i < B->n; i++) {
//     v.x = u.x + B->dx[i];
//     v.y = u.y + B->dy[i];
//     v.z = u.z + B->dz[i];
//     if (ValidVoxel(img,v)){
//       q = GetVoxelIndex(img,v);
//       img->val[q]=YCbCr.val[0];
//       img->Cb[q]=(ushort) YCbCr.val[1];
//       img->Cr[q]=(ushort) YCbCr.val[2];
//     }
//   }

// }

// void DrawLine(Image *img, Voxel u1, Voxel u2, Color YCbCr, AdjRel *B)
// {
//   Point u;
//   Voxel v;
//   int p, k, n;
//   float dx=0, dy=0, dz=0, Dx=u2.x-u1.x, Dy=u2.y-u1.y, Dz=u2.z-u1.z;

//   if ( (!ValidVoxel(img,u1)) || (!ValidVoxel(img,u2)) )
//     Error("Line has end point(s) outside the image domain","DrawLine");

//   if (img->Cb == NULL){
//     img->Cb = AllocUShortArray(img->n);
//     img->Cr = AllocUShortArray(img->n);
//     for (p=0; p < img->n; p++) {
//       img->Cb[p] = 128;
//       img->Cr[p] = 128;
//     }
//   }

//   /* DDA - Digital Differential Analyzer */

//   if (VoxelsAreEqual(u1,u2)) {
//     n = 1;
//   }else{ /* draw line from u1 to u2 */
//     if ((fabs(Dx) >= fabs(Dy))&&(fabs(Dx) >= fabs(Dz))) { /* Dx is the maximum projection of
//                   vector u1u2 */
//       n  = (int)(fabs(Dx)+1);
//       dx = SIGN(Dx);
//       dy = (float)dx*Dy/(float)Dx;
//       dz = (float)dx*Dz/(float)Dx;
//     }else{
//       if ((fabs(Dy) >= fabs(Dx))&&(fabs(Dy) >= fabs(Dz))) { /* Dy is the maximum projection of
//                     vector u1u2 */
//  n  = (int)(fabs(Dy)+1);
//  dy = SIGN(Dy);
//  dx = (float)dy*Dx/(float)Dy;
//  dz = (float)dy*Dz/(float)Dy;
//       } else { /* Dz is the maximum projection of vector u1u2 */
//  n  = (int)(fabs(Dz)+1);
//  dz = SIGN(Dz);
//  dx = (float)dz*Dx/(float)Dz;
//  dy = (float)dz*Dy/(float)Dz;
//       }
//     }
//   }

//   u.x = u1.x;  u.y = u1.y;  u.z = u1.z;
//   for (k=0; k < n; k++) {
//     v.x = ROUND(u.x); v.y = ROUND(u.y); v.z = ROUND(u.z);
//     DrawPoint(img,v,YCbCr,B);
//     u.x = (u.x + dx);
//     u.y = (u.y + dy);
//     u.z = (u.z + dz);
//   }
// }

// void DrawPoints(Image *img, Set *S, Color YCbCr, AdjRel *B)
// {
//   Set *Saux=NULL;
//   Voxel u;
//   int p;

//   if (img->Cb == NULL){
//     img->Cb = AllocUShortArray(img->n);
//     img->Cr = AllocUShortArray(img->n);
//     for (p=0; p < img->n; p++) {
//       img->Cb[p] = 128;
//       img->Cr[p] = 128;
//     }
//   }

//   Saux = S;
//   while (Saux!=NULL) {
//     p   = Saux->elem;
//     u.x = GetXCoord(img,p);
//     u.y = GetYCoord(img,p);
//     u.z = GetZCoord(img,p);
//     DrawPoint(img,u,YCbCr,B);
//     Saux = Saux->next;
//   }

// }

// void DrawBorders(Image *img, Image *label, AdjRel *A, Color YCbCr, AdjRel *B)
// {
//   Voxel u,v;
//   int i,p,q;

//   if ((img->xsize != label->xsize)||
//       (img->ysize != label->ysize)||
//       (img->zsize != label->zsize))
//     Error("Images must have the same domain","DrawBorders");

//   if (img->Cb==NULL)
//     SetCbCr(img,128);

//   for (p=0; p < img->n; p++) {
//     u.x = GetXCoord(label,p);
//     u.y = GetYCoord(label,p);
//     u.z = GetZCoord(label,p);
//     for (i=0; i < A->n; i++) {
//       v.x = u.x + A->dx[i];
//       v.y = u.y + A->dy[i];
//       v.z = u.z + A->dz[i];
//       if (ValidVoxel(label,v)){
//  q = GetVoxelIndex(label,v);
//  if (label->val[p] < label->val[q]){
//    DrawPoint(img,u,YCbCr,B);
//    break;
//  }
//       }
//     }
//   }

// }

// void DrawObject(Image *img, Image *bin, Color YCbCr, AdjRel *B)
// {
//   int p;
//   Voxel u;

//   if ((img->xsize != bin->xsize)||
//       (img->ysize != bin->ysize)||
//       (img->zsize != bin->zsize))
//     Error("Images must have the same domain","DrawObject");

//   SetCbCr(img,128);

//   for (p=0; p < bin->n; p++) {
//     if (bin->val[p]!=0){
//       u.x = GetXCoord(bin,p);
//       u.y = GetYCoord(bin,p);
//       u.z = GetZCoord(bin,p);
//       DrawPoint(img,u,YCbCr,B);
//     }
//   }

// }

// Image *DrawVoxelSamples(DataSet *Z, uchar opt, char *ref_data_type)
// {
//   Image  *img=NULL;
//   int        s,p;
//   Voxel   u;
//   AdjRel *B;
//   Color   RGB, YCbCr;
//   ColorTable *ctb=NULL;

//   if ((Z->ref_data == NULL)||
//       ((strcmp(ref_data_type,"Image")!=0)&&
//        (strcmp(ref_data_type,"MImage")!=0)&&
//        (strcmp(ref_data_type,"FImage")!=0))){
//     Error("Reference data must be an image","DrawVoxelSamples");
//   }

//   if (strcmp(ref_data_type,"Image")==0){
//     Image *aux = (Image *)Z->ref_data;
//     img = CreateImage(aux->xsize,aux->ysize,aux->zsize);
//     SetCbCr(img,128);
//     if (Is3DImage(aux)){
//       MaximumValue(aux);
//       for (p=0; p < aux->n; p++)
//  img->val[p] = (int)(255.*log(aux->val[p]+2)/log(aux->maxval+2));
//     }else{
//       for (p=0; p < aux->n; p++)
//  img->val[p] = aux->val[p];
//     }
//   }else{
//     if (strcmp(ref_data_type,"FImage")==0){
//       FImage *aux = (FImage *)Z->ref_data;
//       img = FImageToImage(aux,255);
//       SetCbCr(img,128);
//     }else{
//       if (strcmp(ref_data_type,"MImage")==0){
//  MImage *aux = (MImage *)Z->ref_data;
//  img = MImageToImage(aux,255,0);
//  SetCbCr(img,128);
//       }
//     }
//   }

//   if (Is3DImage(img))
//     B = Spheric(3.0);
//   else
//     B = Circular(3.0);

//   switch (opt) {

//   case LABEL:

//     ctb = CreateColorTable(MAX(Z->nlabels+2,1000));

//     for (s=0; s < Z->nsamples; s++) {
//       p = Z->sample[s].id;
//       u = GetVoxelCoord(img,p);
//       RGB.val[0] = ctb->color[Z->sample[s].label].val[0];
//       RGB.val[1] = ctb->color[Z->sample[s].label].val[1];
//       RGB.val[2] = ctb->color[Z->sample[s].label].val[2];
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//     }
//     DestroyColorTable(&ctb);

//     break;

//   case CLASS:

//     ctb = CreateColorTable(MAX(Z->nclasses+2,1000));

//     for (s=0; s < Z->nsamples; s++) {
//       p = Z->sample[s].id;
//       u = GetVoxelCoord(img,p);
//       RGB.val[0] = ctb->color[Z->sample[s].class].val[0];
//       RGB.val[1] = ctb->color[Z->sample[s].class].val[1];
//       RGB.val[2] = ctb->color[Z->sample[s].class].val[2];
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//     }
//     DestroyColorTable(&ctb);

//     break;

//   case POINT:

//       RGB.val[0] = 255;
//       RGB.val[1] = 255;
//       RGB.val[2] = 255;
//       YCbCr = RGBtoYCbCr(RGB);
//       for (s=0; s < Z->nsamples; s++) {
//  p = Z->sample[s].id;
//  u = GetVoxelCoord(img,p);
//  DrawPoint(img,u,YCbCr,B);
//       }

//       break;

//     default:

//       Error("Invalid option","DrawVoxelSamples");

//   }

//   DestroyAdjRel(&B);

//   return(img);
// }


// Image *Draw2DFeatureSpace(DataSet *Z, uchar opt, uchar status)
// {
//   Image *img=NULL;
//   int       s,p;
//   Voxel  u;
//   float     xmin,xmax,ymin,ymax,minw,maxw;
//   AdjRel *B=Circular(sqrtf(2.0));
//   Color   RGB, YCbCr;
//   ColorTable *ctb=NULL;
//   float *x=NULL, *y=NULL;

//   if (Z->nfeats != 2)
//     Error("Feature space must be 2D","Draw2DFeatureSpace");

//   minw=xmin=ymin=INFINITY_FLT; maxw=xmax=ymax=-INFINITY_FLT;
//   for (s=0; s < Z->nsamples; s++){
//     if (Z->sample[s].feat[0] < xmin)
//       xmin = Z->sample[s].feat[0];
//     if (Z->sample[s].feat[1] < ymin)
//       ymin = Z->sample[s].feat[1];
//     if (Z->sample[s].feat[0] > xmax)
//       xmax = Z->sample[s].feat[0];
//     if (Z->sample[s].feat[1] > ymax)
//       ymax = Z->sample[s].feat[1];
//     if (Z->sample[s].weight < minw)
//       minw = Z->sample[s].weight;
//     if (Z->sample[s].weight > maxw)
//       maxw = Z->sample[s].weight;
//   }

//   x = AllocFloatArray(Z->nsamples);
//   y = AllocFloatArray(Z->nsamples);

//   for (s=0; s < Z->nsamples; s++){
//     x[s] = 520*(Z->sample[s].feat[0]-xmin)/(xmax - xmin);
//     y[s] = 520*(Z->sample[s].feat[1]-ymin)/(ymax - ymin);
//   }

//   img=CreateImage(540,540,1);
//   SetCbCr(img,128);

//   for (p=0; p < img->n; p++) {
//     RGB.val[0] = 10;
//     RGB.val[1] = 10;
//     RGB.val[2] = 10;
//     YCbCr = RGBtoYCbCr(RGB);
//     img->val[p] = YCbCr.val[0];
//     img->Cb[p]  = YCbCr.val[1];
//     img->Cr[p]  = YCbCr.val[2];
//   }

//   ctb = CreateColorTable(MAX(Z->nclasses+2,1000));

//   for (s=0; s < Z->nsamples; s++) {
//     u.x = (int)x[s] + 5;
//     u.y = (int)y[s] + 5;
//     u.z = 0;
//     p = GetVoxelIndex(img,u);

//     switch (opt) {

//     case WEIGHT:
//       RGB.val[0] = (int) MIN(255.0*(1.0 - (expf(-(Z->sample[s].weight-minw)/(maxw-minw)))),255);
//       RGB.val[1] = 255 - RGB.val[0];
//       RGB.val[2] = 0;
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//       break;

//     case LABEL:
//       RGB.val[0] = ctb->color[Z->sample[s].label].val[0];
//       RGB.val[1] = ctb->color[Z->sample[s].label].val[1];
//       RGB.val[2] = ctb->color[Z->sample[s].label].val[2];
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//       break;

//     case CLASS:

//       RGB.val[0] = ctb->color[Z->sample[s].class].val[0];
//       RGB.val[1] = ctb->color[Z->sample[s].class].val[1];
//       RGB.val[2] = ctb->color[Z->sample[s].class].val[2];
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//       break;

//     case STATUS:
//       if (Z->sample[s].status == status) {
//  RGB.val[0] = ctb->color[Z->sample[s].class].val[0];
//  RGB.val[1] = ctb->color[Z->sample[s].class].val[1];
//  RGB.val[2] = ctb->color[Z->sample[s].class].val[2];
//  YCbCr = RGBtoYCbCr(RGB);
//  DrawPoint(img,u,YCbCr,B);
//       }
//       break;

//     case POINT:
//       RGB.val[0] = 255;
//       RGB.val[1] = 255;
//       RGB.val[2] = 255;
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//       break;

//     default:
//       RGB.val[0] = 255;
//       RGB.val[1] = 255;
//       RGB.val[2] = 255;
//       YCbCr = RGBtoYCbCr(RGB);
//       DrawPoint(img,u,YCbCr,B);
//     }
//   }

//   free(x);
//   free(y);
//   DestroyAdjRel(&B);
//   DestroyColorTable(&ctb);

//   return(img);
// }



// Image *ProjectMeanValue(Image *img, Plane *cutplane, int xviewsize, int yviewsize, int slabsize)
// {
//   Image *mip,*aux;
//   Point  P1,P0;
//   Voxel  u;
//   int       i,xyviewsize,offset,p;

//   /* Initialize output image */

//   mip   = CreateImage(xviewsize,yviewsize,1);

//   /* Starting at slabsize slices before, reslice the scene until
//      slabsize slices after the center of the cut plane and along
//      the direction of v. */

//   xyviewsize = xviewsize*yviewsize;
//   P0   = cutplane->pos;
//   for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
//     P1.x = P0.x + i*cutplane->normal.x;
//     P1.y = P0.y + i*cutplane->normal.y;
//     P1.z = P0.z + i*cutplane->normal.z;
//     u.x  = ROUND(P1.x);
//     u.y  = ROUND(P1.y);
//     u.z  = ROUND(P1.z);

//     if (ValidVoxel(img,u)){
//       SetPlanePos(cutplane,P1.x,P1.y,P1.z);
//       aux  = GetSlice(img,cutplane,xviewsize,yviewsize);
//       for (p=0; p < mip->n; p++) {
//  mip->val[p] += aux->val[p];
//       }
//       DestroyImage(&aux);
//     }
//   }

//   for (p=0; p < mip->n; p++) {
//     mip->val[p] /= (2*slabsize+1);
//   }

//   return(mip);
// }

// Image *ProjectMinValue(Image *img, Plane *cutplane, int xviewsize, int yviewsize, int slabsize)
// {
//   Image *mip,*aux;
//   Point  P1,P0;
//   Voxel  u;
//   int       i,xyviewsize,offset,p;

//   /* Initialize output image */

//   mip   = CreateImage(xviewsize,yviewsize,1);
//   SetImage(mip,MaximumValue(img)+1);

//   /* Starting at slabsize slices before, reslice the scene until
//      slabsize slices after the center of the cut plane and along the
//      direction of v. */

//   xyviewsize = xviewsize*yviewsize;
//   P0 = cutplane->pos;
//   for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
//     P1.x = P0.x + i*cutplane->normal.x;
//     P1.y = P0.y + i*cutplane->normal.y;
//     P1.z = P0.z + i*cutplane->normal.z;
//     u.x  = ROUND(P1.x);
//     u.y  = ROUND(P1.y);
//     u.z  = ROUND(P1.z);

//     if (ValidVoxel(img,u)){
//       SetPlanePos(cutplane,P1.x,P1.y,P1.z);
//       aux  = GetSlice(img,cutplane,xviewsize,yviewsize);
//       for (p=0; p < mip->n; p++) {
//  if (mip->val[p] > aux->val[p])
//    mip->val[p] = aux->val[p];
//       }
//       DestroyImage(&aux);
//     }
//   }

//   for (p=0; p < mip->n; p++) {
//     if (mip->val[p]==(img->maxval+1))
//       mip->val[p]=0;
//   }

//   return(mip);
// }

// Image *ProjectMaxValue(Image *img, Plane *cutplane, int xviewsize, int yviewsize, int slabsize)
// {
//   Image *mip,*aux;
//   Point  P1,P0;
//   Voxel  u;
//   int       i,xyviewsize,offset,p;

//   /* Initialize output image */

//   mip   = CreateImage(xviewsize,yviewsize,1);

//   /* Starting at slabsize slices before, reslice the scene until
//      slabsize slices after the center of the cut plane and along
//      the direction of v. */

//   xyviewsize = xviewsize*yviewsize;
//   P0   = cutplane->pos;
//   for (i=-slabsize,offset=0; i < slabsize; i++,offset+=xyviewsize) {
//     P1.x = P0.x + i*cutplane->normal.x;
//     P1.y = P0.y + i*cutplane->normal.y;
//     P1.z = P0.z + i*cutplane->normal.z;
//     u.x  = ROUND(P1.x);
//     u.y  = ROUND(P1.y);
//     u.z  = ROUND(P1.z);

//     if (ValidVoxel(img,u)){
//       SetPlanePos(cutplane,P1.x,P1.y,P1.z);
//       aux  = GetSlice(img,cutplane,xviewsize,yviewsize);
//       for (p=0; p < mip->n; p++) {
//  if (mip->val[p] < aux->val[p])
//    mip->val[p] = aux->val[p];
//       }
//       DestroyImage(&aux);
//     }
//   }

//   return(mip);
// }

// Image *ProjectMeanObjectValue(Image *img, Image *obj, Plane *cutplane, int xviewsize, int yviewsize, int slabsize)
// {
//   Image  *mip,*aux1,*aux2;
//   Point  P1,P0;
//   Voxel  u;
//   int       i,xyviewsize,offset,p;

//   /* Initialize output image */

//   mip   = CreateImage(xviewsize,yviewsize,1);

//   /* Starting at slabsize slices before, reslice the scene until
//      slabsize slices after the center of the cut plane and along
//      the direction of v. */

//   xyviewsize = xviewsize*yviewsize;
//   P0 = cutplane->pos;
//   for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
//     P1.x = P0.x + i*cutplane->normal.x;
//     P1.y = P0.y + i*cutplane->normal.y;
//     P1.z = P0.z + i*cutplane->normal.z;
//     u.x  = ROUND(P1.x);
//     u.y  = ROUND(P1.y);
//     u.z  = ROUND(P1.z);

//     if (ValidVoxel(img,u)){
//       SetPlanePos(cutplane,P1.x,P1.y,P1.z);
//       aux1  = GetSlice(img,cutplane,xviewsize,yviewsize);
//       aux2  = GetSlice(obj,cutplane,xviewsize,yviewsize);
//  for (p=0; p < mip->n; p++) {
//    if (aux2->val[p] != 0){
//      mip->val[p] += aux1->val[p];
//    }
//  }
//  DestroyImage(&aux1);
//  DestroyImage(&aux2);
//       }
//     }

//   for (p=0; p < mip->n; p++) {
//     mip->val[p]/=(2*slabsize+1);
//   }

//   return(mip);
// }


// Image *ProjectMinObjectValue(Image *img, Image *obj, Plane *cutplane, int xviewsize, int yviewsize, int slabsize)
// {
//   Image  *mip,*aux1,*aux2;
//   Point   P1,P0;
//   Voxel   u;
//   int        i,xyviewsize,offset,p;

//   /* Initialize output image */

//   mip   = CreateImage(xviewsize,yviewsize,1);
//   SetImage(mip,MaximumValue(img)+1);

//   /* Starting at slabsize slices before, reslice the scene until
//      slabsize slices after the center of the cut plane and along
//      the direction of v. */

//   xyviewsize = xviewsize*yviewsize;
//   P0   = cutplane->pos;
//   for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
//     P1.x = P0.x + i*cutplane->normal.x;
//     P1.y = P0.y + i*cutplane->normal.y;
//     P1.z = P0.z + i*cutplane->normal.z;
//     u.x  = ROUND(P1.x);
//     u.y  = ROUND(P1.y);
//     u.z  = ROUND(P1.z);

//     if (ValidVoxel(img,u)){
//       SetPlanePos(cutplane,P1.x,P1.y,P1.z);
//       aux1  = GetSlice(img,cutplane,xviewsize,yviewsize);
//       aux2  = GetSlice(obj,cutplane,xviewsize,yviewsize);
//       for (p=0; p < mip->n; p++) {
//  if (aux2->val[p] != 0){
//    if (mip->val[p] > aux1->val[p])
//      mip->val[p] = aux1->val[p];
//  }
//       }
//       DestroyImage(&aux1);
//       DestroyImage(&aux2);
//     }
//   }

//   for (p=0; p < mip->n; p++) {
//     if (mip->val[p]==(img->maxval+1))
//       mip->val[p]=0;
//   }

//   return(mip);
// }


// Image *ProjectMaxObjectValue(Image *img, Image *obj, Plane *cutplane, int xviewsize, int yviewsize, int slabsize)
// {
//   Image  *mip,*aux1,*aux2;
//   Point  P1,P0;
//   Voxel  u;
//   int       i,xyviewsize,offset,p;

//   /* Initialize output image */

//   mip   = CreateImage(xviewsize,yviewsize,1);

//   /* Starting at slabsize slices before, reslice the scene until
//      slabsize slices after the center of the cut plane and along
//      the direction of v. */

//   xyviewsize = xviewsize*yviewsize;
//   P0   = cutplane->pos;
//   for (i=-slabsize,offset=0; i <= slabsize; i++,offset+=xyviewsize) {
//     P1.x = P0.x + i*cutplane->normal.x;
//     P1.y = P0.y + i*cutplane->normal.y;
//     P1.z = P0.z + i*cutplane->normal.z;
//     u.x  = ROUND(P1.x);
//     u.y  = ROUND(P1.y);
//     u.z  = ROUND(P1.z);

//     if (ValidVoxel(img,u)){
//       SetPlanePos(cutplane,P1.x,P1.y,P1.z);
//       aux1  = GetSlice(img,cutplane,xviewsize,yviewsize);
//       aux2  = GetSlice(obj,cutplane,xviewsize,yviewsize);
//       for (p=0; p < mip->n; p++) {
//  if (aux2->val[p] != 0){
//    if (mip->val[p] < aux1->val[p])
//      mip->val[p] = aux1->val[p];
//  }
//       }
//       DestroyImage(&aux1);
//       DestroyImage(&aux2);
//     }
//   }

//   return(mip);
// }

// Image *CurvilinearSplatting(Image *img, FImage *dist, Plane *cutplane, int viewsize, float depth)
// {
//   Matrix  *R = RotationMatrixToAlignVectorWithZ(cutplane->normal);
//   Image *proj;
//   FImage *zbuff;
//   Voxel u,v;
//   Point P1,P2;
//   int p,q;
//   float dmin=MAX(depth-0.5,0.5),dmax=depth+1.0;


//   /* Initialize output image */

//   proj   = CreateImage(viewsize,viewsize,1);

//   /* Perform voxel splatting, using the direction of the viewing
//      plane's normal and depth information*/

//   zbuff   = CreateFImage(viewsize,viewsize,1);
//   FSetImage(zbuff,INFINITY_FLT);

//   for (u.z=0; u.z < img->zsize; u.z++)
//     for (u.y=0; u.y < img->ysize; u.y++)
//       for (u.x=0; u.x < img->xsize; u.x++){

//  p   = GetVoxelIndex(img,u);

//  if ((dist->val[p] >= dmin)&&(dist->val[p] <= dmax)) {

//    /* Set voxel coordinate to the rotated system with respect
//       to the center of the object */

//    P1.x = u.x - cutplane->pos.x;
//    P1.y = u.y - cutplane->pos.y;
//    P1.z = u.z - cutplane->pos.z;
//    P2   = TransformPoint(R,P1);
//    P1.x = P2.x + viewsize/2.0;
//    P1.y = P2.y + viewsize/2.0;
//    P1.z = P2.z + viewsize/2.0;

//    /* Perform voxel splatting */

//      v.z = 0;

//      v.y = (int) (P1.y);
//      v.x = (int) (P1.x);

//      if (ValidVoxel(proj,v)){
//        q   = GetVoxelIndex(proj,v);
//        if (zbuff->val[q]>P1.z){
//      proj->val[q]  = img->val[p];
//      zbuff->val[q] = P1.z;
//        }
//      }

//      v.y = (int) (P1.y);
//      v.x = (int) (P1.x+1.0);

//      if (ValidVoxel(proj,v)){
//        q   = GetVoxelIndex(proj,v);
//        if (zbuff->val[q]>P1.z){
//      proj->val[q]  = img->val[p];
//      zbuff->val[q] = P1.z;
//        }
//      }

//      v.y = (int) (P1.y+1.0);
//      v.x = (int) (P1.x);

//      if (ValidVoxel(proj,v)){
//        q   = GetVoxelIndex(proj,v);
//        if (zbuff->val[q]>P1.z){
//      proj->val[q]  = img->val[p];
//      zbuff->val[q] = P1.z;
//        }
//      }

//      v.y = (int) (P1.y+1.0);
//      v.x = (int) (P1.x+1.0);

//      if (ValidVoxel(proj,v)){
//        q   = GetVoxelIndex(proj,v);
//        if (zbuff->val[q]>P1.z){
//      proj->val[q]  = img->val[p];
//      zbuff->val[q] = P1.z;
//        }
//      }
//  }
//       }

//   DestroyMatrix(&R);
//   DestroyFImage(&zbuff);

//   return(proj);
// }

// Image *ColorizeComp(Image *label)
// {
//   ColorTable *ctb;
//   Image *img;
//   int p;
//   Color YCbCr;

//   ctb = CreateColorTable(MaximumValue(label)+1);
//   img = CreateImage(label->xsize,label->ysize,label->zsize);
//   SetImage(img,255);
//   SetCbCr(img,128);

//   for (p=0; p < label->n; p++)
//     if (label->val[p]>0){
//       YCbCr=RGBtoYCbCr(ctb->color[label->val[p]]);
//       img->val[p] = YCbCr.val[0];
//       img->Cb[p]  = YCbCr.val[1];
//       img->Cr[p]  = YCbCr.val[2];
//     }
//   DestroyColorTable(&ctb);
//   return(img);
// }

// Image *ColorizeCompOverImage(Image *orig, Image *label)
// {
//   ColorTable *ctb;
//   Image *img;
//   int p;
//   Color YCbCr;

//   if ((orig->xsize != label->xsize)||
//       (orig->ysize != label->ysize)||
//       (orig->zsize != label->zsize))
//     Error("Images must have the same domain","ColorizeCompOnImage");

//   ctb = CreateColorTable(MaximumValue(label)+1);
//   img = CreateImage(label->xsize,label->ysize,label->zsize);
//   for (p=0; p < img->n; p++)
//     img->val[p]=orig->val[p];

//   SetCbCr(img,128);

//   for (p=0; p < label->n; p++)
//     if (label->val[p]>0){
//       YCbCr=RGBtoYCbCr(ctb->color[label->val[p]]);
//       img->val[p] = YCbCr.val[0];
//       img->Cb[p]  = YCbCr.val[1];
//       img->Cr[p]  = YCbCr.val[2];
//     }
//   DestroyColorTable(&ctb);
//   return(img);
// }


