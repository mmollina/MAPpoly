/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2020 Marcelo Mollinari
 
 This file is part of MAPpoly.
 
 MAPpoly is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 For a copy of the GNU General Public License, please visit
 <http://www.gnu.org/licenses/>.
 */

/*
 File: read_mappoly_vcf.cpp

 Description: Set of functions to be used with software R

 Functions used to allow .vcf input file support.

 Functions written by Gabriel Gesteira.

 First version:       2019
 Last update: Sep 23, 2020
 */

#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <zlib.h>
#include <sstream>
#include <string>
using namespace std;
using namespace Rcpp;

//using namespace Rcpp;
// https://www.lemoda.net/c/gzfile-read/

/* Size of the block of memory to use for reading. */
#define LENGTH 4000 // hexadecimel for 16384 or 16.384 KB.

// Getting dosages
int get_dosage(std::string mystring, int gt_pos){
  // mystring is a string containing vcf genotypes
  // gt_pos is the GT position on vcf file
  //  Rcpp::Rcout << "In strsplit" << std::endl;
  char split = ':';
  int flag;
  std::vector<std::string> vec_o_strings;
  int start = 0;
  unsigned int i=0;
  // Looping through string
  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split){
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings.push_back(temp);
      start = i+1;
      i = i+1;
    }
  }
  // Handle the last element
  std::string temp = mystring.substr(start, i - start);
  vec_o_strings.push_back(temp);
  // Testing for correct separation
  std::string test = vec_o_strings[gt_pos-1];
  if(test.size() > 2){
    mystring = vec_o_strings[gt_pos-1];
  }
  start = 0;
  // Checking for NA
  if(mystring[0] == '.'){
    return -1;
  }
  // Checking other formats
  if(mystring[0] == '0' || mystring[0] == '1'){flag = 0;} else {flag = 1;}
  if(flag == 1){
    // Rcpp::Rcout << "Warning: GT field does not have the expected format, please check your file. Returning NA instead.\n";
    return -1;
  }
  // Accounting allele dosages
  for (i = 0; i < mystring.size(); i++){
    if(mystring[i] == '1'){
      start += 1;
    }
  }
  return start;
}

// Getting ploidy information
int get_ploidy(std::string mystring, int gt_pos){
  // mystring is a string containing vcf genotypes
  // gt_pos is the GT position on vcf file
  //  Rcpp::Rcout << "In strsplit" << std::endl;
  char split = ':';
  char split2 = '/';
  char split3 = '|';
  int flag;
  std::vector<std::string> vec_o_strings;
  int start = 0;
  unsigned int i=0;
  // Looping through string
  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split){
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings.push_back(temp);
      start = i+1;
      i = i+1;
    }
  }
  // Handle the last element
  std::string temp = mystring.substr(start, i - start);
  vec_o_strings.push_back(temp);
  // Testing for correct separation
  std::string test = vec_o_strings[gt_pos-1];
  if(test.size() > 2){
    mystring = vec_o_strings[gt_pos-1];
  }
  start = 0;
  // // Checking for NA (not needed for ploidy)
  // if(mystring[0] == '.'){
  //   return -1;
  // }
  // Checking other formats
  if(mystring[0] == '0' || mystring[0] == '1' || mystring[0] == '.'){flag = 0;} else {flag = 1;}
  if(flag == 1){
    // Rcpp::Rcout << "Warning: GT field does not have the expected format, please check your file. Returning NA instead.\n";
    return -1;
  }
  // Looping through string2
  start = 0;
  for(i = 0; i < mystring.size(); i++){
    if( mystring[i] != split2 && mystring[i] != split3){
      start += 1;
    }
  }  
  return start;
}

// Getting depth information
int get_depth(std::string mystring, int dp_pos){
  // mystring is a string containing vcf information for each marker-ind
  // dp_pos is the DP position on vcf file
  //  Rcpp::Rcout << "In strsplit" << std::endl;
  char split = ':';
  std::vector<std::string> vec_o_strings;
  int start = 0;
  unsigned int i=0;
  // Looping through string
  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split){
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings.push_back(temp);
      start = i+1;
      i = i+1;
    }
  }
  // Handle the last element
  std::string temp = mystring.substr(start, i - start);
  vec_o_strings.push_back(temp);
  // Returning DP
  if (vec_o_strings.size() < static_cast<std::vector<std::string>::size_type>(dp_pos)){
    start = 0;
  } else if (vec_o_strings[dp_pos-1] == ".") {
    start = 0;
  } else {
    start = stoi(vec_o_strings[dp_pos-1]);
  }
  return start;
}

// Getting PL information
std::vector<double> get_probabilities(std::string mystring, int pl_pos){
  // mystring is a string containing vcf genotypes
  // gt_pos is the GT position on vcf file
  //  Rcpp::Rcout << "In strsplit" << std::endl;
  char split = ':';
  char split2 = ',';
  std::vector<std::string> vec_o_strings;
  std::vector<std::string> vec_o_strings2;
  std::vector<double> final_vec_out2(1);
  int start = 0, size = 0, s = 0;
  double sum = 0.0;
  unsigned int i=0;
  // Looping through string
  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split){
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings.push_back(temp);
      start = i+1;
      i = i+1;
    }
  }
  // Handle the last element
  std::string temp = mystring.substr(start, i - start);
  vec_o_strings.push_back(temp);
  // Checking PL information
  if (vec_o_strings.size() < static_cast<std::vector<std::string>::size_type>(pl_pos)){
    final_vec_out2[0] = -1.0;
    return final_vec_out2;
  }
  // Testing for correct separation
  std::string test = vec_o_strings[pl_pos-1];
  if(test.size() > 2){
    mystring = vec_o_strings[pl_pos-1];
  } else {
    mystring = "0,0";
  }
  // Looping through selected PL string
  start = 0;
  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split2){
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings2.push_back(temp);
      start = i+1;
      i = i+1;
    }
  }
  // Handle the last element
  temp = mystring.substr(start, i - start);
  vec_o_strings2.push_back(temp);
  size = vec_o_strings2.size();
  std::vector<double> final_vec(size);
  std::vector<double> final_vec_out(size);
  // Converting to double and transforming scale (from phred to probabilities)
  for (s=0; s < size; s++){
    final_vec[s] = stod(vec_o_strings2[s]);
    final_vec[s] = final_vec[s]/(-10);
    final_vec[s] = exp(final_vec[s]);
    sum += final_vec[s];
  }
  // Rescale PLs to sum=1
  if (sum > 0.001){
    for (s=0; s < size; s++){
      final_vec_out[s] = final_vec[s]/sum;
    }
  } else {
    for (s=0; s < size; s++){
      final_vec_out[s] = final_vec[s];
    }
  }
  sum = 0.0;
  return final_vec_out;
}

// Getting PL for all markers and all individuals
// [[Rcpp::export(name=".vcf_get_probabilities")]]
Rcpp::List vcf_get_probabilities(Rcpp::StringMatrix& mat, int pl_pos){
  int rowmat = mat.nrow();
  int colmat = mat.ncol();
  int k, l, s, size;
  Rcpp::List results(rowmat);
  Rcpp::NumericMatrix partial_results;
  std::vector<double> first, out;
  std::vector<double> get_probabilities(std::string mystring, int pl_pos);
  // Looping through matrix
  for (k=0; k < rowmat; k++){
    // Rcout << "The value of k : " << k << "\n";
    first = get_probabilities(as<std::string>(mat(k,0)), pl_pos);
    size = first.size();
    Rcpp::NumericMatrix partial_results(colmat,size);
    for (l=0; l < colmat; l++){
      // Rcout << "The value of l : " << l << "\n";
      out = get_probabilities(as<std::string>(mat(k,l)), pl_pos);
      if (out.size() == static_cast<std::vector<double>::size_type>(size)){
        for (s=0; s < size; s++){
          // Rcout << "The value of s : " << s << "\n";
          partial_results(l,s) = out[s];
        }
      } else {
        for (s=0; s < size; s++){
          partial_results(l,s) = -1.0;
        }
      }
      s=0;
    }
    results(k) = partial_results;
  }
  return results;
}

// Getting dosages for all markers and all individuals
// [[Rcpp::export(name=".vcf_transform_dosage")]]
Rcpp::NumericMatrix vcf_transform_dosage(Rcpp::StringMatrix& mat, int gt_pos){
  int rowmat = mat.nrow();
  int colmat = mat.ncol();
  int k, l;
  Rcpp::NumericMatrix results(rowmat,colmat);
  int get_dosage(std::string mystring, int gt_pos);
  for (k=0; k < rowmat; k++){
    for (l=0; l < colmat; l++){
      results(k,l) = get_dosage(as<std::string>(mat(k,l)), gt_pos);
    }
  }
  return results;
}

// Getting ploidy for all markers and all individuals
// [[Rcpp::export(name=".vcf_get_ploidy")]]
Rcpp::NumericMatrix vcf_get_ploidy(Rcpp::StringMatrix& mat, int gt_pos){
  int rowmat = mat.nrow();
  int colmat = mat.ncol();
  int k, l;
  Rcpp::NumericMatrix results(rowmat,colmat);
  int get_ploidy(std::string mystring, int gt_pos);
  for (k=0; k < rowmat; k++){
    for (l=0; l < colmat; l++){
      results(k,l) = get_ploidy(as<std::string>(mat(k,l)), gt_pos);
    }
  }
  return results;
}

// Getting depth for all markers and all individuals
// [[Rcpp::export(name=".vcf_get_depth")]]
Rcpp::NumericMatrix vcf_get_depth(Rcpp::StringMatrix& mat, int dp_pos){
  int rowmat = mat.nrow();
  int colmat = mat.ncol();
  int k, l;
  Rcpp::NumericMatrix results(rowmat,colmat);
  int get_dosage(std::string mystring, int gt_pos);
  for (k=0; k < rowmat; k++){
    for (l=0; l < colmat; l++){
      results(k,l) = get_depth(as<std::string>(mat(k,l)), dp_pos);
    }
  }
  return results;
}

// Rcpp::NumericMatrix vcf_transform_dosage(Rcpp::StringMatrix& mat, int gt_pos){
//   int rowmat = mat.nrow();
//   int colmat = mat.ncol();
//   int k, l;
//   int res;
//   char split = ':';
//   int flag;
//   std::vector<std::string> vec_o_strings;
//   int start = 0;
//   int i=0;
//   Rcpp::NumericMatrix results(rowmat,colmat);
//   std::string mystring;
//   
//   for (k=0; k < rowmat; k++){
//     for (l=0; l < colmat; l++){
//       start = 0;
//       // Defining mystring
//       mystring = mat(k,l);
//       // Rcpp::Rcout << mystring << '\n';
//       // Looping through string
//       for(i = 1; i < mystring.size(); i++){
//         if( mystring[i] == split){
//           std::string temp = mystring.substr(start, i - start);
//           vec_o_strings.push_back(temp);
//           start = i+1;
//           i = i+1;
//         }
//       }
//       
//       // Handle the last element
//       std::string temp = mystring.substr(start, i - start);
//       vec_o_strings.push_back(temp);
//       
//       // Testing for correct separation
//       std::string test = vec_o_strings[gt_pos-1];
//       if(test.size() > 2){
//         mystring = vec_o_strings[gt_pos-1];
//       }
//       start = 0;
//       
//       // Checking other formats
//       if(mystring[0] == '0' || mystring[0] == '1'){flag = 0;} else {flag = 1;}
//       
//       // Checking for NA
//       if(mystring[0] == '.'){
//         start = -1;
//       } else if(flag == 0){
//         // Accounting allele dosages
//         for (i = 0; i < mystring.size(); i++){
//           if(mystring[i] == '1'){
//             start += 1;
//           }
//         }
//       }
//       else if(flag == 1){
//         Rcpp::Rcout << "Warning: GT field does not have the expected format, please check your file. Returning NA instead.\n";
//         start = -1;
//       }
// 
//       results(k,l) = start;
//     }
//   }
//   return results;
// }
