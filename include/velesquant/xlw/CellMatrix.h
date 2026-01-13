#pragma once
#include <vector>
#include <iostream>
#include <string>

namespace velesquant {

namespace xlw {

    class CellMatrix {
    public:
        struct Cell {
            double val;
            Cell(double v = 0.0) : val(v) {}
            double NumericValue() const { return val; }
            operator double() const { return val; }
            bool IsANumber() const { return true; }
            Cell& operator+=(double v) { val += v; return *this; }
            Cell& operator=(double v) { val = v; return *this; }
        };

        CellMatrix() {}
        CellMatrix(size_t rows, size_t cols) {
            data_.resize(rows, std::vector<Cell>(cols));
        }
        CellMatrix(double v) {
            data_.resize(1, std::vector<Cell>(1, Cell(v)));
        }
        CellMatrix(const std::vector<double>& v) {
            data_.resize(v.size(), std::vector<Cell>(1));
            for(size_t i=0; i<v.size(); ++i) data_[i][0] = Cell(v[i]);
        }
        CellMatrix(const std::string& s) {
             // Basic support for string if passed, though Cell is double-only in my stub.
             // If legacy code passes strings, this might be issue.
             // For now assume double-only or use 0.0.
              data_.resize(1, std::vector<Cell>(1, 0.0));
        }

        Cell& operator()(size_t i, size_t j) { 
            if (i >= data_.size() || j >= data_[0].size()) {
                 if (i >= data_.size()) data_.resize(i+1);
                 if (j >= data_[i].size()) data_[i].resize(j+1);
            }
            return data_[i][j]; 
        }
        const Cell& operator()(size_t i, size_t j) const { return data_[i][j]; }
        
        size_t RowsInStructure() const { return data_.size(); }
        size_t ColumnsInStructure() const { return data_.empty() ? 0 : data_[0].size(); }
        
        // Aliases for compatibility
        size_t rows() const { return RowsInStructure(); }
        size_t columns() const { return ColumnsInStructure(); }

    private:
        std::vector<std::vector<Cell>> data_;
    };
    
    typedef CellMatrix MyMatrix;
    typedef std::vector<double> MyArray;

}
}