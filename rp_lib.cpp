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
 *
 * rp_lib.cpp : modified version of rp.cpp for library-use (originally written by Nicola Prezza)
 *
 *  Created on: Sep 3, 2025
 *      Author: Hiroki Shibata
 *
 */

#include "re_pair.hpp"

#include <chrono>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <stack>
#include <fstream>
#include <cmath>
#include <cassert>

#include <lf_queue.hpp>
#include <ll_vec.hpp>
#include <ll_el.hpp>

#include "internal/hf_queue.hpp"
#include "internal/skippable_text.hpp"
#include "internal/text_positions.hpp"
#include "internal/packed_gamma_file3.hpp"

namespace rp {
  namespace detail {

    template<bool Use64> struct Types;

    template<> struct Types<false> {
      using text_t = skippable_text32_t;
      using TP_t   = text_positions32_t;
      using hf_q_t = hf_queue32_t;
      using lf_q_t = lf_queue32_t;
      using itype  = uint32_t;
    };

    template<> struct Types<true> {
      using text_t = skippable_text64_t;
      using TP_t   = text_positions64_t;
      using hf_q_t = hf_queue64_t;
      using lf_q_t = lf_queue64_t;
      using itype  = uint64_t;
    };

    template<bool Use64>
    struct CodecImpl {
      using T      = Types<Use64>;
      using text_t = typename T::text_t;
      using TP_t   = typename T::TP_t;
      using hf_q_t = typename T::hf_q_t;
      using lf_q_t = typename T::lf_q_t;
      using itype  = typename T::itype;
      using cpair  = typename hf_q_t::cpair;

      itype X=0;
      itype last_freq = 0;
      itype n_distinct_freqs = 0;

      std::vector<itype> A;
      std::vector<std::pair<itype, itype> > G;
      std::vector<itype> T_vec;

      void new_high_frequency_queue(hf_q_t & Q, TP_t & TP, text_t & T, uint64_t min_freq){
        itype j = 0;
        itype n = TP.size();
        int old_perc = 0;
        int perc;
        itype n_pairs = 0;
        while(j<n){
          itype k = 1;
          while( j<TP.size()-1 &&
                 T.pair_starting_at(TP[j]) != T.blank_pair() &&
                 T.pair_starting_at(TP[j]) == T.pair_starting_at(TP[j+1]) ){
            j++;
            k++;
          }
          if(k>=min_freq){
            n_pairs++;
          }
          j++;
        }
        itype max_d = 256+T.size()/min_freq;
        Q.init(max_d,min_freq);
        j = 0;
        while(j<n){
          itype P_ab = j;
          itype k = 1;
          cpair ab = T.blank_pair();
          while( j<TP.size()-1 &&
                 T.pair_starting_at(TP[j]) != T.blank_pair() &&
                 T.pair_starting_at(TP[j]) == T.pair_starting_at(TP[j+1]) ){
            ab = T.pair_starting_at(TP[j]);
            j++;
            k++;
          }
          if(k>=min_freq){
            Q.insert({ab, P_ab, k, k});
          }
          j++;
        }
      }

      template<typename queue_t>
      void synchronize(queue_t & Q, TP_t & TP, text_t & T, cpair AB){
        assert(Q.contains(AB));
        auto q_el = Q[AB];
        itype P_AB = q_el.P_ab;
        itype L_AB = q_el.L_ab;
        itype F_AB = q_el.F_ab;

        itype freq_AB = 0;

        assert(P_AB+L_AB <= TP.size());
        TP.cluster(P_AB,P_AB+L_AB);
        assert(TP.is_clustered(P_AB,P_AB+L_AB));

        itype j = P_AB;
        while(j<P_AB+L_AB){
          itype p = j;
          itype k = 1;
          cpair XY = T.pair_starting_at(TP[j]);
          while( j<(P_AB+L_AB)-1 &&
                 XY != T.blank_pair() &&
                 XY == T.pair_starting_at(TP[j+1]) ){
            j++;
            k++;
          }
          freq_AB = XY == AB ? k : freq_AB;
          if(k >= Q.minimum_frequency()){
            if(XY != AB){
              assert(XY != T.blank_pair());
              assert(not Q.contains(XY));
              Q.insert({XY,p,k,k});
              assert(TP.contains_only(p,p+k,XY));
            }else if(XY == AB){
              assert(Q.contains(AB));
              Q.update({AB,p,k,k});
              assert(TP.contains_only(p,p+k,AB));
            }
          }
          j++;
        }

        assert(Q.contains(AB));
        if(freq_AB < Q.minimum_frequency()){
          Q.remove(AB);
        }
        assert(not Q.contains(AB) || Q[AB].F_ab == Q[AB].L_ab);
      }

      template<typename queue_t>
      void synchro_or_remove_pair(queue_t & Q, TP_t & TP, text_t & T, cpair ab){
        assert(Q.contains(ab));
        auto q_el = Q[ab];
        itype F_ab = q_el.F_ab;
        itype L_ab = q_el.L_ab;

        if(F_ab <= L_ab/2){
          synchronize<queue_t>(Q, TP, T, ab);
        }else{
          if(F_ab < Q.minimum_frequency()){
            Q.remove(ab);
          }
        }
      }

      uint64_t wd(uint64_t x){
        auto w = 64 - __builtin_clzll(uint64_t(x));
        return x == 0 ? 1 : w;
      }

      const std::pair<itype,itype> nullpair_v = {~itype(0),~itype(0)};

      template<typename queue_t>
      uint64_t substitution_round(queue_t & Q, TP_t & TP, text_t & T){
        using ctype = typename text_t::char_type;
        cpair AB = Q.max();
        G.push_back(AB);
        assert(Q.contains(AB));
        assert(Q[AB].F_ab >= Q.minimum_frequency());
        auto q_el = Q[AB];
        itype F_AB = q_el.F_ab;
        itype P_AB = q_el.P_ab;
        itype L_AB = q_el.L_ab;

        uint64_t f_replaced = F_AB;

        n_distinct_freqs += (F_AB != last_freq);
        last_freq = F_AB;

        for(itype j = P_AB; j<P_AB+L_AB;++j){
          itype i = TP[j];
          if(T.pair_starting_at(i) == AB){
            ctype A = AB.first;
            ctype B = AB.second;
            cpair xA = T.pair_ending_at(i);
            cpair By = T.next_pair(i);
            assert(xA == T.blank_pair() or xA.second == A);
            assert(By == T.blank_pair() or By.first == B);
            T.replace(i,X);
            assert(By == T.blank_pair() || T.pair_starting_at(i) == cpair(X,By.second));
            if(Q.contains(xA) && xA != AB){
              Q.decrease(xA);
            }
            if(Q.contains(By) && By != AB){
              Q.decrease(By);
            }
          }
        }

        for(itype j = P_AB; j<P_AB+L_AB;++j){
          itype i = TP[j];
          assert(T.pair_starting_at(i) != AB);
          if(T[i] == X){
            cpair xX = T.pair_ending_at(i);
            cpair Xy = T.pair_starting_at(i);
            ctype A = AB.first;
            ctype B = AB.second;
            ctype x = xX.first == X ? B : xX.first;
            ctype y = Xy.second == X ? A : Xy.second;
            cpair xA = xX == T.blank_pair() ? xX : cpair {x,A};
            cpair By = Xy == T.blank_pair() ? Xy : cpair {B,y};
            if(Q.contains(By) && By != AB){
              synchro_or_remove_pair<queue_t>(Q, TP, T, By);
            }
            if(Q.contains(xA) && xA != AB){
              synchro_or_remove_pair<queue_t>(Q, TP, T, xA);
            }
          }
        }

        assert(Q.contains(AB));
        synchronize<queue_t>(Q, TP, T, AB);
        assert(not Q.contains(AB));
        X++;
        return f_replaced;
      }

      void compute_repair(const std::vector<uint8_t>& in_bytes){
        double alpha = 0.66;
        uint64_t B = 50;

        itype n = in_bytes.size();
        itype sigma = 0;

        itype min_high_frequency = 0;
        itype lf_queue_capacity = 0;

        const itype null = ~itype(0);

        std::vector<itype> char_to_int(256,null);

        min_high_frequency = std::pow(  n, alpha  );
        min_high_frequency = min_high_frequency <2 ? 2 : min_high_frequency;

        // std::cout << "File size = " << n << " characters"  << std::endl;
        // std::cout << "cut-off frequency = " << min_high_frequency  << std::endl;

        itype width = 64 - __builtin_clzll(uint64_t(n));
        itype max_d = 256+n/min_high_frequency;

        text_t T(n);

        itype j = 0;

        // std::cout << "filling skippable text with text characters ... " << std::flush;

        for(auto c : in_bytes){
          if(char_to_int[uint8_t(c)] == null){
            char_to_int[uint8_t(c)] = sigma;
            A.push_back(uint8_t(c));
            sigma++;
          }
          T.set(j++,char_to_int[uint8_t(c)]);
        }

        // std::cout << "done. " << std::endl << std::endl;
        // std::cout << "alphabet size is " << sigma  << std::endl << std::endl;
        // std::cout << "initializing and sorting text positions vector ... " << std::flush;

        TP_t TP(&T,min_high_frequency);

        // std::cout << "done. Number of text positions containing a high-frequency pair: " << TP.size() << std::endl;

        X =  sigma;

        // std::cout << "\nSTEP 1. HIGH FREQUENCY PAIRS" << std::endl << std::endl;
        // std::cout << "inserting pairs in high-frequency queue ... " << std::flush;

        hf_q_t HFQ;
        new_high_frequency_queue(HFQ, TP, T, min_high_frequency);

        // std::cout << "done. Number of distinct high-frequency pairs = " << HFQ.size() << std::endl;
        // std::cout << "Replacing high-frequency pairs ... " << std::endl;

        int last_perc = -1;
        uint64_t F = 0;

        while(HFQ.max() != HFQ.nullpair()){
          auto f = substitution_round<hf_q_t>(HFQ, TP, T);
          if(last_perc == -1){
            F = f;
            last_perc = 0;
          }else{
            int perc = 100-(100*f)/F;
            if(perc > last_perc+4){
              last_perc = perc;
              // std::cout << perc << "%" << std::endl;
            }
          }
        }

        // std::cout << "done. " << std::endl;
        // std::cout << "Peak queue size = " << HFQ.peak() << " (" << double(HFQ.peak())/double(n) << "n)" << std::endl;

        // std::cout << "\nSTEP 2. LOW FREQUENCY PAIRS" << std::endl << std::endl;
        // std::cout << "Re-computing TP array ... " << std::flush;

        TP.fill_with_text_positions();

        // std::cout << "done." << std::endl;
        // std::cout << "Sorting  TP array ... " << std::flush;
        TP.cluster();
        // std::cout << "done." << std::endl;

        // std::cout << "Counting low-frequency pairs ... " << std::flush;

        uint64_t n_lf_pairs = 0;

        uint64_t f = 1;
        for(uint64_t i=1;i<TP.size();++i){
          if(T.pair_starting_at(TP[i]) == T.pair_starting_at(TP[i-1])){
            f++;
          }else{
            f=1;
            n_lf_pairs++;
          }
        }
        // std::cout << "done. Number of distict low-frequency pairs: "<< n_lf_pairs << std::endl;

        // std::cout << "Filling low-frequency queue ... " << std::flush;

        lf_q_t LFQ(min_high_frequency-1);

        f = 1;
        using el_t = typename lf_q_t::el_type;

        for(uint64_t i=1;i<TP.size();++i){
          if(T.pair_starting_at(TP[i]) == T.pair_starting_at(TP[i-1])){
            f++;
          }else{
            if(f>1){
              cpair ab = T.pair_starting_at(TP[i-1]);
              assert(i>=f);
              itype P_ab = i - f;
              itype L_ab = f;
              itype F_ab = f;
              el_t el = {ab,P_ab,L_ab,F_ab};
              LFQ.insert(el);
            }
            f=1;
          }
        }

        // std::cout << "done." << std::endl;
        // std::cout << "Replacing low-frequency pairs ... " << std::endl;

        std::pair<itype,itype> replaced = {0,0};

        last_perc = -1;
        uint64_t tl = T.number_of_non_blank_characters();

        while(LFQ.max() != LFQ.nullpair()){
          auto f2 = substitution_round<lf_q_t>(LFQ, TP, T);
          int perc = 100-(100*T.number_of_non_blank_characters())/tl;
          if(perc>last_perc+4){
            last_perc = perc;
            // std::cout << perc << "%" << std::endl;
          }
        }

        // std::cout << "done. " << std::endl;
        // std::cout << "Peak queue size = " << LFQ.peak() << " (" << double(LFQ.peak())/double(n) << "n)" << std::endl;
        // std::cout << "Compressing grammar and storing it to file ... " << std::endl << std::endl;

        for(itype i=0;i<T.size();++i){
          if(not T.is_blank(i)) T_vec.push_back(T[i]);
        }
      }

      RPRePair compute_repair(std::string in){
        std::ifstream ifs(in);
        char c;
        std::vector<uint8_t> in_bytes;
        while(ifs.get(c)){
          in_bytes.emplace_back(c);
        }
        compute_repair(in_bytes);
        return {A, G, T_vec};
      }

      void decompress_vec(std::vector<itype> & A_, std::vector<std::pair<itype,itype> > & G_, std::vector<itype> & Tc_, std::ofstream & ofs){
        std::stack<itype> S;
        std::string buffer;
        int buf_size = 1000000;

        for(itype i = 0;i<Tc_.size();++i){
          S.push(Tc_[i]);
          while(!S.empty()){
            itype X_ = S.top();
            S.pop();
            if(X_<A_.size()){
              char c = A_[X_];
              buffer.push_back(c);
              if(buffer.size()==buf_size){
                ofs.write(buffer.c_str(),buffer.size());
                buffer = std::string();
              }
            }else{
              auto ab = G_[X_-A_.size()];
              S.push(ab.second);
              S.push(ab.first);
            }
          }
        }
        if(buffer.size()>0) ofs.write(buffer.c_str(),buffer.size());
      }
    };

  } // namespace detail

  RPRePair compute_repair(const std::vector<uint8_t>& in_bytes) {
    detail::CodecImpl<false> impl;
    impl.compute_repair(in_bytes);
    return {
      impl.A, impl.G, impl.T_vec
    };
  }

  void compress_and_write(const std::vector<uint8_t>& in_bytes, const std::string& out){
    detail::CodecImpl<false> impl;
    // std::cout << "Compressing file " << in << std::endl;
    // std::cout << "Output will be saved to file " << out << std::endl << std::endl;
    impl.compute_repair(in_bytes);
    packed_gamma_file3<typename detail::CodecImpl<false>::itype> out_file(out);
    out_file.compress_and_store(impl.A, impl.G, impl.T_vec);
  }

  void compress_file(const std::string& in, const std::string& out) {
    detail::CodecImpl<false> impl;
    // std::cout << "Compressing file " << in << std::endl;
    // std::cout << "Output will be saved to file " << out << std::endl << std::endl;
    impl.compute_repair(in);
    packed_gamma_file3<typename detail::CodecImpl<false>::itype> out_file(out);
    out_file.compress_and_store(impl.A, impl.G, impl.T_vec);
  }

  void decompress_file(const std::string& in, const std::string& out) {
    // std::cout << "Decompressing archive " << in << std::endl;
    // std::cout << "Output will be saved to " << out << std::endl;
    auto pgf = packed_gamma_file3<>(in, false);

    {
      detail::CodecImpl<false> impl32;
      std::vector<typename detail::CodecImpl<false>::itype> A;
      std::vector<std::pair<typename detail::CodecImpl<false>::itype,typename detail::CodecImpl<false>::itype> > G;
      std::vector<typename detail::CodecImpl<false>::itype> Tc;
      pgf.read_and_decompress(A,G,Tc);
      std::ofstream ofs(out, std::ios::binary);
      impl32.decompress_vec(A,G,Tc,ofs);
      ofs.close();
    }
    // std::cout << "done." << std::endl;
  }

} // namespace re_pair
