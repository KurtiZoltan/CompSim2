#include "Utils.hpp"
#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <cmath>

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

template <typename T>
T f(T x, T y)
{
    return std::exp(-(x * x + y * y) / 10);
}

template <typename T>
T gauss(T x, T y)
{
    T sigma2 = 5;
    return std::exp(-(x * x + y * y) / (2 * sigma2));
}

template <typename T>
T dgauss(T x, T y)
{
    T sigma2 = 5;
    return 1 / sigma2 * ((x * x + y * y) / sigma2 - 2) * std::exp(-(x * x + y * y) / (2 * sigma2));
}

template <typename T>
T g(T x, T y)
{
    return -5 * std::sin(2 * x - y);
}

template <typename T>
void fillGrid(T (*f)(T x, T y), T* grid, u32 n, u32 m)
{
    assert(n == m & "test function works on a square");
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            T x = static_cast<T>(i) / (n - 1);
            T y = static_cast<T>(j) / (m - 1);
            grid[j * n + i] = f(x, y);
        }
    }
}

template <typename T>
void laplace(T* gridIn, T* gridOut, u32 n, u32 m, T h, T* kernel, u32 kernelRadius)
{
    //bulk calculations
    for (u32 i = kernelRadius; i < n - kernelRadius; i++)
    {
        for (u32 j = kernelRadius; j < m - kernelRadius; j++)
        {
            T temp = 0;
            for (i32 dx = -static_cast<i32>(kernelRadius); dx <= static_cast<i32>(kernelRadius); dx++)
            {
                for (i32 dy = -static_cast<i32>(kernelRadius); dy <= static_cast<i32>(kernelRadius); dy++)
                {
                    temp += kernel[(dy + kernelRadius) * (2 * kernelRadius + 1) + (dx + kernelRadius)] * gridIn[(j + dy) * n + i + dx];
                }
            }
            gridOut[j * n + i] = temp / h / h;
        }
    }
}

template <typename T>
void printGrid(T* grid, u32 n, u32 m, const std::string fileName)
{
    std::fstream file;
    file.open(fileName, std::ios::out);
    if (!file.is_open())
    {
        std::cout << "File '" << fileName << "' couldn't be opened.\n";
        return;
    }
    file << std::setprecision(16);
    for (u32 i = 0; i < n; i++)
    {
        for (u32 j = 0; j < m; j++)
        {
            file << grid[j * n + i] << " ";
        }
        file << "\n";
    }
}

template <typename T>
T calcError(T* grid1, T* grid2, u32 N)
{
    T error = 0;
    for (u32 i = 2; i < N - 2; i++)
    {
        for (u32 j = 2; j < N - 2; j++)
        {
            error += std::abs(grid1[j * N + i] - grid2[j * N + i]);
        }
    }
    error /= (N - 4) * (N - 4);
    return error;
}

template <typename T>
void test(u32 N, std::fstream& file)
{
    T h = 1.0 / (N - 1);
    T* testFunction = allocateGrid<T>(N, N);
    T* exact = allocateGrid<T>(N, N);
    T* numeric = allocateGrid<T>(N, N);
    fillGrid(gauss, testFunction, N, N);
    fillGrid(dgauss, exact, N, N);
    
    T kernel1a[] = {
        1.0 / 6,  2.0 / 3, 1.0 / 6,
        2.0 / 3,-10.0 / 3, 2.0 / 3,
        1.0 / 6,  2.0 / 3, 1.0 / 6,
    };
    T kernel1b[] = {
        0,  1, 0,
        1, -4, 1,
        0,  1, 0,
    };
    u32 kernel1Radius = 1;
    
    T kernel2a[] = {
         0.0/60, -2.0/60,   -1.0/60, -2.0/60,  0.0/60,
        -2.0/60, 16.0/60,   52.0/60, 16.0/60, -2.0/60,
        -1.0/60, 52.0/60, -252.0/60, 52.0/60, -1.0/60,
        -2.0/60, 16.0/60,   52.0/60, 16.0/60, -2.0/60,
         0.0/60, -2.0/60,   -1.0/60, -2.0/60,  0.0/60,
    };
    u32 kernel2Radius = 2;
    
    laplace(testFunction, numeric, N, N, h, kernel1a, kernel1Radius);
    T errora = calcError(numeric, exact, N);
    
    laplace(testFunction, numeric, N, N, h, kernel1b, kernel1Radius);
    T errorb = calcError(numeric, exact, N);
    
    laplace(testFunction, numeric, N, N, h, kernel2a, kernel2Radius);
    T errorc = calcError(numeric, exact, N);
    
    file << N << " " << std::setprecision(16) << errora << " " << errorb << " " << errorc << "\n";
    
    freeGrid(numeric);
    freeGrid(exact);
    freeGrid(testFunction);
}

typedef f32 floatType;

int main()
{
    std::fstream file64;
    std::string name64 = "./2dlaplaceDouble.txt";
    file64.open(name64, std::ios::out);
    if (!file64.is_open())
    {
        std::cout << name64 << " couldn't be opened.\n";
        return 1;
    }
    std::fstream file32;
    std::string name32 = "./2dlaplaceFloat.txt";
    file32.open(name32, std::ios::out);
    if (!file32.is_open())
    {
        std::cout << name32 << " couldn't be opened.\n";
        return 1;
    }
    
    u32 Nmin = 10;
    u32 Nmax = 10000;
    u32 numPoints = 50;
    for (u32 t = 0; t < numPoints; t++)
    {
        u32 N = Nmin * std::exp(std::log(static_cast<f64>(Nmax) / Nmin) * t / (numPoints - 1));
        std::cout << N << "\n";
        test<f64>(N, file64);
        test<f32>(N, file32);
    }
    return 0;
}
