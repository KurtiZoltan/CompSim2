#include "Utils.hpp"
#include <iostream>
#include <cmath>
#include <utility>

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

void initRegion(u8* region, bool (*f)(u32 i, u32 j), u32 n, u32 m)
{
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            region[j * m + i] = f(i, j) ? 1 : 0;
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
    u8 type; // 0 - Dirichlet, 1 - Neumann
    u32 i, j;
    union 
    {
        f64 normal[2];
        f64 value;
    };
};

std::pair<BoundaryElement*, u32> allocateBoundary(u8* region, u32 n, u32 m, BoundaryElement (*boundary)(u32 i, u32 j), f64 dx, f64 dy)
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
                boundaryElements[boundaryCount] = boundary(i, j);
                boundaryCount++;
            }
        }
    }
    return std::make_pair(boundaryElements, boundaryCount);
}

void freeBoundary(BoundaryElement* normals)
{
    delete[] normals;
}

void laplace35(f64* gridIn, u8* region, f64* gridOut, u32 n, u32 m, f64 h)
{
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            f64 temp = 0;
            if (region[j * n + i] == 2)
            {
                temp = (
                    (gridIn[(j - 1) * n + i - 1] + gridIn[(j - 1) * n + i + 1] + gridIn[(j + 1) * n + i - 1] + gridIn[(j + 1) * n + i + 1]) + 
                    (gridIn[j * n + i - 1] + gridIn[j * n + i + 1] + gridIn[(j - 1) * n + i] + gridIn[(j + 1) * n + i]) * 4 +
                    -(gridIn[j * n + i]) * 20
                    ) / 6;
            }
            else if (region[j * n + i] > 2)
            {
                temp = (
                    -(gridIn[j * n + i - 2] + gridIn[j * n + i + 2] + gridIn[(j - 2) * n + i] + gridIn[(j + 2) * n + i]) * 1 +
                    -(gridIn[(j + 1) * n + i - 2] + gridIn[(j + 1) * n + i + 2] + gridIn[(j - 2) * n + i + 1] + gridIn[(j + 2) * n + i + 1] +
                      gridIn[(j - 1) * n + i - 2] + gridIn[(j - 1) * n + i + 2] + gridIn[(j - 2) * n + i - 1] + gridIn[(j + 2) * n + i - 1]) * 2 +
                    (gridIn[(j - 1) * n + i - 1] + gridIn[(j - 1) * n + i + 1] + gridIn[(j + 1) * n + i - 1] + gridIn[(j + 1) * n + i + 1]) * 16 + 
                    (gridIn[j * n + i - 1] + gridIn[j * n + i + 1] + gridIn[(j - 1) * n + i] + gridIn[(j + 1) * n + i]) * 52 +
                    -(gridIn[j * n + i]) * 252
                    ) / 60;
            }
            gridOut[j * n + i] = temp / h / h;
        }
    }
}

void laplace3(f64* gridIn, u8* region, f64* gridOut, u32 n, u32 m, f64 h)
{
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            f64 temp = 0;
            if (region[j * n + i] > 1)
            {
                temp = (
                    (gridIn[(j - 1) * n + i - 1] + gridIn[(j - 1) * n + i + 1] + gridIn[(j + 1) * n + i - 1] + gridIn[(j + 1) * n + i + 1]) + 
                    (gridIn[j * n + i - 1] + gridIn[j * n + i + 1] + gridIn[(j - 1) * n + i] + gridIn[(j + 1) * n + i]) * 4 +
                    -(gridIn[j * n + i]) * 20
                    ) / 6;
            }
            gridOut[j * n + i] = temp / h / h;
        }
    }
}

void laplace3ani(f64* gridIn, u8* region, f64* gridOut, u32 n, u32 m, f64 h)
{
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            f64 temp = 0;
            if (region[j * n + i] > 1)
            {
                temp = 
                    (gridIn[j * n + i - 1] + gridIn[j * n + i + 1] + gridIn[(j - 1) * n + i] + gridIn[(j + 1) * n + i]) +
                    -(gridIn[j * n + i]) * 4;
            }
            gridOut[j * n + i] = temp / h / h;
        }
    }
}

void applyBoundary(f64* grid, BoundaryElement* boundary, u32 boundaryCount)
{
    
}

//-------------------------------------------------------------------------------------------------

u32 N = 1024;

bool circle(u32 i, u32 j)
{
    f64 x = 2 * (static_cast<f64>(i) / (N - 1) - 0.5);
    f64 y = 2 * (static_cast<f64>(j) / (N - 1) - 0.5);
    f64 R = 0.9;
    return x * x + y * y < R * R;
}

BoundaryElement circleBoundary(u32 i, u32 j)
{
    f64 nx = 2 * (static_cast<f64>(i) / (N - 1) - 0.5);
    f64 ny = 2 * (static_cast<f64>(j) / (N - 1) - 0.5);
    f64 norm = std::sqrt(nx * nx + ny * ny);
    nx /= norm;
    ny /= norm;
    BoundaryElement ret;
    ret.type = 1;
    ret.i = i;
    ret.j = j;
    ret.normal[0] = nx;
    ret.normal[1] = ny;
    return ret;
}

int main()
{
    f64 dx = 1.0 / (N - 1);
    f64* T = allocateGrid<f64>(N, N);
    f64* laplace = allocateGrid<f64>(N, N);
    u8* region = allocateGrid<u8>(N, N);
    initRegion(region, circle, N, N);
    auto ret = allocateBoundary(region, N, N, circleBoundary, dx, dx);
    BoundaryElement* boundary = ret.first;
    u32 boundaryCount = ret.second;
    
    
    
    freeBoundary(boundary);
    freeGrid(region);
    freeGrid(laplace);
    freeGrid(T);
    return 0;
}