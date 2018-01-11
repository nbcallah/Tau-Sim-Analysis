/*  Header file for Sebastiano Vigna's xorshift1024*φ
    I'm keeping things as close to the original as possible*/

#include <stdint.h>

#ifndef GRANDPARENT_H
#define GRANDPARENT_H

uint64_t s[16];
int p;

uint64_t next(void);
void jump(void);
double nextU01();
void initxorshift();

#endif