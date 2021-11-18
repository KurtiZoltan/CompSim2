#include "Utils.hpp"
#include <iostream>
#include <cmath>

//-------------------------------------------------------------------------------------------------
//helper functions
template <typename T>
T* allocateGrid(u32 n, u32 m)
{
    T* grid = new T[n * m];
    for (u32 i = 0; i < n * m; i++)
    {
        grid[i] = 0;
    }
    return grid;
}

template <typename T>
void freeGrid(T* grid)
{
    delete[] grid;
}

void fillRegion(u8* region, bool (*f)(u32 i, u32 j, u32 n, u32 m), u32 n, u32 m)
{
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            region[j * m + i] = f(i, j, n, m) ? 1 : 0;
        }
    }
    for (u32 d = 2; 2 * d < n && 2 * d < m && d < 256; d++)
    {
        for (u32 i = d - 1; i < n - d + 1; i++)
        {
            for (u32 j = d - 1; j < m - d + 1; j++)
            {
                bool interior = true;
                for (i32 di = -1; di <= 1; di++)
                {
                    for (i32 dj = -1; dj <= 1; dj++)
                    {
                        if (region[(j + dj) * n + i + di] < d - 1)
                            interior = false;
                    }
                }
                if (interior)
                    region[j * n + i] = d;
            }
        }
    }
}

struct BoundaryElement
{
    u32 i, j;
    f64 nx, ny;
};

BoundaryElement* allocateBoundaryNormals(u8* region, u32 n, u32 m, void (*normal)(u32 i, u32 j, u32 n, u32 m, f64 dx, f64 dy, f64& nx, f64& ny), f64 dx, f64 dy)
{
    u32 boundaryCount = 0;
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            if (region[j * n + i] == 1)
            {
                boundaryCount++;
            }
        }
    }
    std::cout << boundaryCount << " boundary elements.\n";
    BoundaryElement* boundaryElements = new BoundaryElement[boundaryCount];
    boundaryCount = 0;
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            if (region[j * n + i] == 1)
            {
                boundaryElements[boundaryCount].i = i;
                boundaryElements[boundaryCount].j = j;
                normal(i, j, n, m, dx, dy, boundaryElements[boundaryCount].nx, boundaryElements[boundaryCount].ny);
                boundaryCount++;
            }
        }
    }
    return boundaryElements;
}

void freeBoundaryNormals(BoundaryElement* normals)
{
    delete[] normals;
}

//-------------------------------------------------------------------------------------------------

bool circle(u32 i, u32 j, u32 n, u32 m)
{
    f64 x = 2 * (static_cast<f64>(i) / (n - 1) - 0.5);
    f64 y = 2 * (static_cast<f64>(j) / (m - 1) - 0.5);
    f64 R = 0.9;
    return x * x + y * y < R * R;
}

void circleNormal(u32 i, u32 j, u32 n, u32 m, f64 dx, f64 dy, f64& nx, f64& ny)
{
    f64 x = 2 * (static_cast<f64>(i) / (n - 1) - 0.5);
    f64 y = 2 * (static_cast<f64>(j) / (m - 1) - 0.5);
    nx = -x / dx;
    ny = -y / dy;
    f64 norm = std::sqrt(nx * nx + ny * ny);
    nx /= norm;
    ny /= norm;
}

int main()
{
    u32 N = 1024;
    f64 h = 1.0 / (N - 1);
    f64* T = allocateGrid<f64>(N, N);
    u8* region = allocateGrid<u8>(N, N);
    fillRegion(region, circle, N, N);
    BoundaryElement* normals = allocateBoundaryNormals(region, N, N, circleNormal, h, h);
    
    
    
    freeBoundaryNormals(normals);
    freeGrid(region);
    freeGrid(T);
    return 0;
}