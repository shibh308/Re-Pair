/*
 *  This file is part of Re-Pair.
 *  Copyright (c) by
 *  Nicola Prezza <nicola.prezza@gmail.com>, Philip Bille, and Inge Li Gørtz
 *
 *   Re-Pair is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   Re-Pair is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 *
 * rp_pair.hpp : header file for Re-Pair library (added by Hiroki Shibata)
 *
 *  Created on: Sep 3, 2025
 *      Author: Hiroki Shibata
 *
 */

#ifndef GRAMMAR_WITH_LZSE_ENCODING_RE_PAIR_HPP
#define GRAMMAR_WITH_LZSE_ENCODING_RE_PAIR_HPP

#include <string>
#include <vector>

namespace rp {
  struct RPRePair {
    using itype = uint32_t;
    std::vector<itype> A;
    std::vector<std::pair<itype, itype> > G;
    std::vector<itype> T_vec;
  };
  RPRePair compute_repair(const std::vector<uint8_t>& in_bytes);
  void encode_repair(RPRePair rp, const std::string& out);
  RPRePair decode_repair(const std::string& in);
  std::vector<uint8_t> reconstruct(const RPRePair& rp);
} // namespace re_pair


#endif //GRAMMAR_WITH_LZSE_ENCODING_RE_PAIR_HPP
