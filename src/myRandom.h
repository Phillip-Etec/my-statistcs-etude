#ifndef READANDWRITE_H_INCLUDED
#define READANDWRITE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <time.h>


unsigned int stringToInt(std::string bits)
{
    int byte;
    byte = stoi(bits, 0, 2);

    return byte;
}

unsigned int randomSeed()
{
    std::string randomByte;
    std::ifstream text ("2020-03-22.txt");
    int randomSeed;
    getline(text, randomByte);
    text.close();
    unsigned int a;
    srand (time(NULL));

    a = rand() % 8388577;

    randomByte = randomByte.substr(a, 31);
    randomSeed = stringToInt(randomByte);

    return randomSeed;
}

int actualRandom(int minValue, int maxValue)
{
    srand(time(NULL));
    int trulyRandom = rand();
    return trulyRandom;
}
//math:sqrt(-2 * math:log(local:random())) * math:sin(2 * math:pi() * local:random())

#endif // READANDWRITE_H_INCLUDED
