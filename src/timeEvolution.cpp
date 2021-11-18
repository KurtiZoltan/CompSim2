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

std::pair<BoundaryElement*, u32> allocateBoundary(u8* region, u32 n, u32 m, BoundaryElement (*boundary)(u32 i, u32 j))
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

void applyBoundary(f64* grid, u8* region, u32 n, u32 m, BoundaryElement* boundary, u32 boundaryCount)
{
    //https://folk.ntnu.no/leifh/teaching/tkt4140/._main056.html
    for (u32 k = 0; k < boundaryCount; k++)
    {
        u32 i = boundary[k].i;
        u32 j = boundary[k].j;
        switch (boundary[k].type)
        {
            case 0:
                grid[j * n + i] = boundary[k].value;
                break;
            case 1:
                f64 nx = boundary[k].normal[0];
                f64 ny = boundary[k].normal[1];
                f64 sum = 0;
                f64 weightSum = 0;
                for (i32 di = -1; di <= 1; di++)
                {
                    for (i32 dj = -1; dj <= 1; dj++)
                    {
                        if (region[(j + dj) * n + i + di] > 1)
                        {
                            f64 weight = 1 / (di * nx + dj * ny) * (di * nx + dj * ny);
                            sum += weight * grid[(j + dj) * n + i + di];
                            weightSum += weight;
                        }
                    }
                }
                if (weightSum < 1e-3)
                {
                    std::cout << "Problem with applying Neumann boundary condition at " << i << " " << j << ".\n";
                    std::cout << "weightSum=" << weightSum << "\n";
                    std::cout << "nx=" << nx << " ny=" << ny << "\n";
                    for (i32 di = -1; di <= 1; di++)
                    {
                        for (i32 dj = -1; dj <= 1; dj++)
                        {
                            std::cout << "di=" << di << " dj=" << dj << "region=" << static_cast<u32>(region[(j + dj) * n + i + di]);
                            if (region[(j + dj) * n + i + di] > 1)
                            {
                                f64 weight = 1.0 / (di * nx + dj * ny) * (di * nx + dj * ny);
                                std::cout << " weight=" << weight << "\n";
                            }
                            std::cout << "\n";
                        }
                    }
                    std::cout << "\n";
                }
                grid[j * n + i] = sum / weightSum;
                break;
        }
    }
}

void printTemperature(f64* T, u32 n, u32 m, f64 min, f64 max, const std::string& name)
{
    //https://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            f64 t = T[j * n + i];
            f64 ratio = 2 * (t - min) / (max - min);
            u8 b = static_cast<u8>(std::max(0.0, 255 * (1 - ratio)));
            u8 r = static_cast<u8>(std::max(0.0, 255 * (ratio - 1)));
            u8 g = 255 - r - b;
            
        }
    }
}
//-------------------------------------------------------------------------------------------------

u32 N = 1024;
f64 R = 0.9;
f64 x_0 = 0;
f64 y_0 = -0.5;
f64 sigma = 0.1;
f64 Tfinal = 10;

bool circle(u32 i, u32 j)
{
    //R = 0.9 circle in an NxN square of points representing the region [-1,1] x [-1,1]
    f64 x = 2 * (static_cast<f64>(i) / (N - 1) - 0.5);
    f64 y = 2 * (static_cast<f64>(j) / (N - 1) - 0.5);
    return x * x + y * y < R * R;
}

BoundaryElement circleBoundary(u32 i, u32 j)
{
    //Neumann type boundary condition for the circle
    f64 nx = 2 * (static_cast<f64>(i) / (N - 1) - 0.5);
    f64 ny = 2 * (static_cast<f64>(j) / (N - 1) - 0.5);
    f64 norm = std::sqrt(nx * nx + ny * ny);
    nx /= norm;
    ny /= norm;
    BoundaryElement ret;
    ret.type = 1; //Neumann type type
    ret.i = i;
    ret.j = j;
    ret.normal[0] = nx;
    ret.normal[1] = ny;
    return ret;
}

int main()
{
    f64 dx = 1.0 / (N - 1);
    f64 dt = dx;
    //allocating arrays
    f64* T = allocateGrid<f64>(N, N); //temperature
    f64* laplace = allocateGrid<f64>(N, N); //storage for the laplace of the temperature
    u8* region = allocateGrid<u8>(N, N); //the region in which to solve the equation
    //init the region based on the circle(i, j) < 0 condition
    initRegion(region, circle, N, N);
    //crete the boundary conditions for the boundary points determined by initRegion
    auto ret = allocateBoundary(region, N, N, circleBoundary);
    BoundaryElement* boundary = ret.first;
    u32 boundaryCount = ret.second;
    //init temperature with a unit gauss (x0, y0, sigma) times Tfinal / Area
    for (u32 i = 0; i < N; i++)
    {
        for (u32 j = 0; j < N; j++)
        {
            if (region[j * N + i] > 0)
            {
                f64 x = 2 * (static_cast<f64>(i) / (N - 1) - 0.5);
                f64 y = 2 * (static_cast<f64>(j) / (N - 1) - 0.5);
                T[j * N + i] = Tfinal / (M_PI * R * R) * 1 / (2 * M_PI * sigma * sigma) * std::exp(- ((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0)) / (2 * sigma * sigma));
            }
        } 
    }
    
    //time evolution
    u32 frame = 0;
    for (f64 t = 0; t < 1; t += dt, frame++)
    {
        std::cout << "Printing frame #" << frame << " time=" << t << "\n";
        printTemperature(T, N, N);
        
        laplace35(T, region, laplace, N, N, dx);
        //laplace3(T, region, laplace, N, N, dx);
        //laplace3ani(T, region, laplace, N, N, dx);
        
        for (u32 i = 0; i < N; i++)
        {
            for (u32 j = 0; j < N; j++)
            {
                T[j * N + i] += laplace[j * N + i] * dt;
            }
        }
        applyBoundary(T, region, N, N, boundary, boundaryCount);
    }
    
    freeBoundary(boundary);
    freeGrid(region);
    freeGrid(laplace);
    freeGrid(T);
    return 0;
}