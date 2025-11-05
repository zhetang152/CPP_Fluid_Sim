#ifndef GRID_H
#define GRID_H

#include<vector>
#include<cassert>
#include<algorithm>
#include<cstddef> //for size_t

template<typename T>    
class Grid{
private:
    int m_width;
    int m_height;
    int m_depth;
    std::vector<T> m_data;
public:
    using sizeType = std::size_t;

    Grid(int width, int height,int depth,const T& initial_value=T()):m_width(width), m_height(height), m_depth(depth){
        //对维度的合法性检查
        assert(width > 0 && "Grid width must be positive");
        assert(height > 0 && "Grid height must be positive");
        assert(depth > 0 && "Grid depth must be positive");

        sizeType totalSize = static_cast<sizeType>(width) * height * depth;
        m_data.resize(totalSize, initial_value);
        }

    sizeType getIndex(int i, int j, int k) const {
        assert(i >= 0 && i < m_width && "Grid access out of bounds (i)");
        assert(j >= 0 && j < m_height && "Grid access out of bounds (j)");
        assert(k >= 0 && k < m_depth && "Grid access out of bounds (k)");

        return static_cast<sizeType>(i) + 
               static_cast<sizeType>(j) * m_width +
               static_cast<sizeType>(k) * m_width * m_height;
    }
    T& operator()(int i, int j, int k) {
        // 使用三维索引公式
        return m_data[getIndex(i, j, k)];
    }
    const T& operator()(int i, int j, int k) const {
        return m_data[getIndex(i, j, k)];
    }

    void fill(const T& value) {
        std::fill(m_data.begin(), m_data.end(), value);
    }
    int getWidth() const noexcept { return m_width; }
    int getHeight() const noexcept { return m_height; }
    int getDepth() const noexcept { return m_depth; }
    
};
#endif // GRID_H