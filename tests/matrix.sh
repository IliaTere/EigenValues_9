#!/bin/bash
# 4 runs
echo $'\n' "a.txt --------------"
./a.out 4 0 1e-13 0 matr/a.txt;
./a.out 4 0 1e-14 0 matr/a.txt;
./a.out 4 0 1e-15 0 matr/a.txt;
echo $'\n' "a20.txt --------------"
./a.out 4 0 1e-13 0 matr/a20.txt
./a.out 4 0 1e-14 0 matr/a20.txt
./a.out 4 0 1e-15 0 matr/a20.txt
echo $'\n' "a40.txt --------------"
./a.out 4 0 1e-13 0 matr/a40.txt
./a.out 4 0 1e-14 0 matr/a40.txt
./a.out 4 0 1e-15 0 matr/a40.txt
echo $'\n' "b.txt --------------"
./a.out 4 0 1e-13 0 matr/b.txt
./a.out 4 0 1e-14 0 matr/b.txt
./a.out 4 0 1e-15 0 matr/b.txt
echo $'\n' "c.txt --------------"
./a.out 6 0 1e-13 0 matr/c.txt
./a.out 6 0 1e-14 0 matr/c.txt
./a.out 6 0 1e-15 0 matr/c.txt
echo $'\n' "d.txt --------------"
./a.out 6 0 1e-13 0 matr/d.txt
./a.out 6 0 1e-14 0 matr/d.txt
./a.out 6 0 1e-15 0 matr/d.txt
echo $'\n' "e.txt --------------"
./a.out 6 0 1e-13 0 matr/e.txt
./a.out 6 0 1e-14 0 matr/e.txt
./a.out 6 0 1e-15 0 matr/e.txt
echo $'\n' "f.txt --------------"
./a.out 4 0 1e-14 0 matr/f.txt
echo $'\n' "g.txt --------------"
./a.out 4 0 1e-14 0 matr/g.txt
