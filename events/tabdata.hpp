// vim: set ai et ts=2 sw=2 tw=80:
//

#pragma once

#include <fmt/format.h>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

using VecD = std::vector<double>;
using VecStr = std::vector<std::string>;
using DataMap = std::map<std::string, VecD>;

class TabData {
public:
  explicit TabData(const std::string &file_name, const VecStr &keys,
                   const bool skip_zero_row = true)
      : _keys{keys} {
    std::iota(_keys_index.begin(), _keys_index.end(), 0);
    std::ifstream input_stream(file_name, std::ios::in);
    while (input_stream.good()) {
      for (auto it = keys.cbegin(); it != keys.cend();) {
        double tmp{0.0};
        input_stream >> tmp;

        // If there is a blank line before EOF, ignore it and break out
        if (input_stream.eof())
          break;

        // ignore entire line if it starts with non-numeric data
        if (it == keys.begin()) {
          if (input_stream.rdstate() == std::ios::failbit) {
            input_stream.clear();
            input_stream.ignore(std::numeric_limits<std::streamsize>::max(),
                                '\n');
            continue;
          }

          if (skip_zero_row) {
            // Skip entire row if first element is 0
            if (std::abs(tmp) < std::numeric_limits<double>::min()) {
              input_stream.ignore(std::numeric_limits<std::streamsize>::max(),
                  '\n');
              continue;
            }
          }
        }
        _data[*it].push_back(tmp);
        // fmt::print("{:15.5E}", _data[*it].back());
        // if (keys.front() == *it) {
        //   fmt::println("{:g}", _data[*it].back());
        // }
        ++it;
      }
      // fmt::println("{:s}", "");
    }
  }

  VecD col_data(const unsigned int column_num) const {
    auto key = _keys.at(column_num);
    return _data.at(key);
  }

  VecD col_data(const std::string &key) const {
    VecD res;
    return _data.at(key);
  }

  inline VecStr keys() const { return _keys; }

  inline DataMap get_data() const { return _data; }

  void print() const {
    for (auto r : std::views::iota(0u, _data.at(_keys.front()).size())) {
      for (auto k : _keys) {
        fmt::print("{:15.5E}", _data.at(k).at(r));
      }
      fmt::print("{:s}", "\n");
    }
  }

private:
  VecStr _keys;
  std::vector<unsigned int> _keys_index;
  DataMap _data;
};
