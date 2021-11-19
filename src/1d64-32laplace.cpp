#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

f64 L = 2 * PI;

template<typename T>
void laplace(const T* state, T* out, u32 N, T h)
{
    T h2 = h * h;
    out[0] = (-2 * state[0] + state[1]) / h2;
    out[N-1] = (-2 * state[N-1] + state[N-2]) / h2;
    for (u32 i = 1; i < N - 1; i++)
    {
        out[i] = (state[i-1] - 2 * state[i] + state[i+1]) / h2;
    }
}

template<typename T>
void advancedLaplace(const T* state, T* out, u32 N, T h)
{
    T h2 = h * h;
    
    out[0] = (-2 * state[0] + state[1]) / h2;
    out[N-1] = (-2 * state[N-1] + state[N-2]) / h2;
    
    out[1] = (state[0] - 2 * state[1] + state[2]) / h2;
    out[N-2] = (state[N-3] - 2 * state[N-2] + state[N-1]) / h2;
    
    for (u32 i = 2; i < N - 2; i++)
    {
        T estimate1 = state[i-1] - 2 * state[i] + state[i+1];
        T estimate2 = (state[i-2] - 2 * state[i] + state[i+2]) / 4;
        out[i] = (4 * estimate1 - estimate2) / 3 / h2;
    }
}

template<typename T>
void test(u32 N, std::ostream& os)
{
    T* state = new T[N];
    T* derivative1 = new T[N];
    T* derivative2 = new T[N];
    T h = L / (N - 1);
    
    for (u32 i = 0; i < N; i++)
    {
        state[i] = std::sin(i * h);
    }
    
    laplace(state, derivative1, N, h);
    advancedLaplace(state, derivative2, N, h);
    
    T max1 = 0;
    T max2 = 0;
    for (u32 i = 2; i < N - 2; i++)
    {
        if (std::abs(state[i] + derivative1[i]) > max1)
        {
            max1 = std::abs(state[i] + derivative1[i]);
        }
        if (std::abs(state[i] + derivative2[i]) > max2)
        {
            max2 = std::abs(state[i] + derivative2[i]);
        }
    }
    
    os << N << " " << max1 << " " << max2 << "\n";
    
    delete[] derivative2;
    delete[] derivative1;
    delete[] state;
}

int main()
{
    u32 pointNum = 1000;
    u32 minPoints = 18;
    u32 maxPoints = 50000;
    f32 scale = std::log(maxPoints / minPoints) / pointNum;
    std::fstream file;
    
    file.open("./laplace32.txt", std::ios::out);
    for (u32 i = 0; i < pointNum; i++)
    {
        test<f32>(static_cast<u32>(minPoints * std::exp(i * scale)), file);
    }
    file.close();
    
    file.open("./laplace64.txt", std::ios::out);
    for (u32 i = 0; i < pointNum; i++)
    {
        test<f64>(static_cast<u32>(minPoints * std::exp(i * scale)), file);
    }
    file.close();
    return 0;
}