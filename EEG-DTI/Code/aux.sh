#!/bin/bash          

mkdir network

cd tmp_data

matlab compute_similarity.m > output.out
