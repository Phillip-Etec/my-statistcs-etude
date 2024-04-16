#ifndef MYSTATS_H_INCLUDED
#define MYSTATS_H_INCLUDED

#include <math.h>
#include <bits/stdc++.h>
#include <random>
#include <utility>

#include "myRandom.h"

int randal(int min, int max) //range : [min, max]
{
   static bool first = true;
   if (first)
   {
      srand( time(NULL) ); //seeding for the first time only!
      first = false;
   }
   return min + rand() % (( max + 1 ) - min);
}

double randfrom(double min, double max)
{
    static bool first = true;

    if(first)
    {
        srand( time(NULL) ); //seeding for the first time only!
        first = false;
    }
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int biggestIndex(int data[], int length)
{
    int i, largestIndex=0, previous = 0;
    int greaterValue = data[0], prevgrtr = data[0];

    for (i = 1; i < length; i++)
    {
        greaterValue = data[i]*(data[i] > greaterValue)+greaterValue*(greaterValue>=data[i]);
        largestIndex = i*(data[i] > prevgrtr)+previous*(prevgrtr >= data[i]);
        previous = largestIndex;
        prevgrtr = greaterValue;
    }

    return largestIndex;
}

int biggestValue(int data[], int length)
{
    int i;
    int greaterValue = data[0];

    for (i = 1; i < length; i++)
    {
        greaterValue = data[i]*(data[i] > greaterValue)+greaterValue*(greaterValue>=data[i]);
    }

    return greaterValue;
}

double biggestValue(double data[], int length)
{
    int i;
    double greaterValue = data[0];

    for (i = 1; i < length; i++)
    {
        greaterValue = data[i]*(data[i] > greaterValue)+greaterValue*(greaterValue>=data[i]);
    }

    return greaterValue;
}

int smallestValue(int data[], int length)
{
    int i;
    int smallerValue = data[0];

    for (i = 1; i < length; i++)
    {
        smallerValue = data[i]*(data[i] < smallerValue)+smallerValue*(smallerValue<=data[i]);
    }

    return smallerValue;
}

double smallestValue(double data[], int length)
{
    int i;
    double smallerValue = data[0];

    for (i = 1; i < length; i++)
    {
        smallerValue = data[i]*(data[i] < smallerValue)+smallerValue*(smallerValue<=data[i]);
    }

    return smallerValue;
}

double minimumAmplitude(double data[], int length)
{
    double smol=smallestValue(data, length), byg = biggestValue(data, length);
    double amplitude = byg - smol;

    return amplitude / 4;
}

double setAmplitude(double data[], int length, int intervals)
{
    double smol=smallestValue(data, length), byg = biggestValue(data, length);
    double amplitude = byg - smol;

    return amplitude / intervals;
}

int intervals(double data[], int length, double amplitude)
{
    double difference = (biggestValue(data, length) - smallestValue(data, length));
    return (difference / amplitude) + 1;
}

int intervals(double data[], int length, double amplitude, bool inclusive)
{
    double difference = (biggestValue(data, length) - smallestValue(data, length));
    return (difference / amplitude) + 1*inclusive;
}

int absoluteFrequency(int data[], int length, int value)
{
    int repeats = 0;

    for (int i = 0; i < length; i++)
    {
        repeats += 1*(data[i]==value);
    }
    return repeats;
}

int absoluteFrequency(double data[], int length, double amplitude, double value)
{
    int repeats = 0;
    double smol = smallestValue(data, length), lowerBound, upperBoound, valuesInterval;
    valuesInterval = (value - smol) / amplitude;
    lowerBound = (int) valuesInterval * amplitude + smol;
    upperBoound = (int) (valuesInterval+1) * amplitude + smol;

    for (int i = 0; i < length; i++)
    {
        repeats += 1*(data[i]>=lowerBound || data[i]<=upperBoound);
    }
    return repeats;
}

double mean(int data[], int length)
{
    int i, sum = 0;

    double mean;

    for (i = 0; i < length; i++)
    {
      sum += data[i];
    }
    mean = (double) sum / length;

    return mean;
}

double mean(double data[], int length) // no need for amplitude support
{
    int i;
    double sum = 0;
    double mean;

    for (i = 0; i < length; i++)
    {
      sum += data[i];
    }
    mean = (double) sum / length;

    return mean;
}

double median(int data[], int length)
{
    std::sort(data, data+length);
    int middle = ((ceil(double(length)/2))-1)+1*(length%2==0), middlow = (floor(double(length)/2))-1;

    return (double(data[middle])+double((data[middlow])*(length%2==0)))/(1+(length%2==0));
}

double median(double data[], int length)
{
    std::sort(data, data+length);
    int middle = ((ceil(double(length)/2))-1)+1*(length%2==0), middlow = (floor(double(length)/2))-1;

    return (double(data[middle])+double((data[middlow])*(length%2==0)))/(1+(length%2==0));
}

int mode(int data[], int length)
{
    bool indexesToMiss[length];
    int indexCount[length];
    int i=0, j;
    for(i=0; i<length; i++)
    {
        indexesToMiss[i] = false;
        indexCount[i] = 1;
    }
    for(i=0; i<length; i++)
    {
        if(indexesToMiss[i])
            continue;
        else
        {
            for(j=i+1*(i<length-1); j<length; j++)
            {
                indexesToMiss[j]=(data[i] == data[j]);
                indexCount[i] += 1*(data[i] == data[j]);
            }
        }
    }
    int biggest = biggestIndex(indexCount, *(&indexCount + 1) - indexCount);
    return data[biggest];
}

double mode(double data[], int length, double amplitude) // support ok I guess too
{
    int intervlao = intervals(data, length, amplitude);
    int indexCount[intervlao];
    int i=0;
    double currentInterval;
    double minim = smallestValue(data, length);
    for(i=0; i<intervlao; i++)
    {
        indexCount[i] = 1;
    }
    for(i=0; i<length; i++)
    {
        currentInterval = (data[i] - minim) / amplitude;

        indexCount[(int)currentInterval] += 1;
    }
    int biggest = biggestIndex(indexCount, *(&indexCount + 1) - indexCount);

    return amplitude * biggest + smallestValue(data, length);
}

int range(int data[], int length)
{
    return biggestValue(data, length) - smallestValue(data, length)+1;
}

double range(double data[], int length)
{
    return biggestValue(data, length) - smallestValue(data, length);
}

double standardDeviation(int data[], int length)
{
    double ma = mean(data, length);
    int i;
    double squaredResults[length];
    for(i=0; i<length; i++)
    {
        squaredResults[i] = pow(double(data[i]-ma), 2);
    }
    ma = mean(squaredResults, length);
    double ans = double(sqrt(ma));
    return ans;
}

double standardDeviation(double data[], int length)
{
    double ma = mean(data, length);
    int i;
    double squaredResults[length];
    for(i=0; i<length; i++)
    {
        squaredResults[i] = pow(double(data[i]-ma), 2);
    }
    ma = mean(squaredResults, length);
    double ans = double(sqrt(ma));
    return ans;
}

double variance(int data[], int length)
{
    double me = mean(data, length);
    int i;
    double squaredDifferences[length];
    for(i=0; i<length; i++)
    {
        squaredDifferences[i] = pow(double(data[i]-me), 2);
    }
    double ans = mean(squaredDifferences, length);
    return ans;
}

double variance(double data[], int length)
{
    double me = mean(data, length);
    int i;
    double squaredDifferences[length];
    for(i=0; i<length; i++)
    {
        squaredDifferences[i] = pow(double(data[i]-me), 2);
    }
    double ans = mean(squaredDifferences, length);
    return ans;
}

double percentileOf(int data[], int length, int percentOf)
{
    int i, quantityBelow = 0, ans;
    for(i=0; i<length; i++)
    {
        quantityBelow += 1*(data[i] < percentOf);
    }
    ans = double(quantityBelow)/length * 100;
    return ans;
}

double percentileOf(double data[], int length, double percentOf)
{
    int i, quantityBelow = 0, ans;
    for(i=0; i<length; i++)
    {
        quantityBelow += 1*(data[i] < percentOf);
    }
    ans = double(quantityBelow)/length * 100;
    return ans;
}

int percentile(int data[], int length, int percent)
{
    int percentileIndex = (int)floor(length/100.0*percent);
    std::sort(data, data+length);
    int thePercentile = data[percentileIndex];

    return thePercentile;
}

int percentile(int data[], int length, double percent)
{
    int percentileIndex = (int)floor(length/100.0*percent);
    std::sort(data, data+length);
    int thePercentile = data[percentileIndex];

    return thePercentile;
}

double percentile(double data[], int length, int percent)
{
    int percentileIndex = (int)floor(length/100.0*percent);
    std::sort(data, data+length);
    int thePercentile = data[percentileIndex];

    return thePercentile;
}

double percentile(double data[], int length, double percent)
{
    int percentileIndex = (int)floor(length/100.0*percent);
    std::sort(data, data+length);
    int thePercentile = data[percentileIndex];

    return thePercentile;
}

void showHistogramBad(int data[], int length)
{
    int i, j, beggining, finish, height;
    beggining = smallestValue(data, length);
    finish = biggestValue(data, length);
    height = finish-beggining;
    int hist[(height+1)];
    int count[(*(&hist +1)-hist)];

    //initialize histogram as every number from lowest of data to biggest of data
    for(i=0, j=beggining; i<(*(&hist +1)-hist); i++, j++)
    {
        hist[i] = j;
    }

    //initialize count as 0 for each
    for(i=0; i<(*(&count + 1)-count); i++)
    {
        count[i] = 0;
    }

    //for every element in count, assign the quantity of times an element of data appears
    for(i = 0; i<(*(&count +1)-count); i++)
    {
        for(j=0; j<length; j++)
        {
            count[i] += 1*(hist[i]==data[j]);
        }
    }

    for(i=0, j=beggining; i<(*(&count +1)-count); i++, j++)
    {
        std::cout << j << '\t';
        for(int ammount=0; ammount<count[i]; ammount++)
        {
            std::cout << '*';
        }
        std::cout << std::endl;
    }
}

int histchar(int m, int n, int k)
{
    int flat = 124, down = 47, up = 98, dive = 60, soar = 62 ;

    int desiredChar = flat*(m==n && k==m) ^
                      down*(k<m && n<k) ^
                        up*((k>m && n>k)) ^
                        dive*(k<m && k<n) ^
                        soar*(k>m && k>m) //^
                        //flat*(k==m || k==n)
                        ;

    return desiredChar;
}

int histchar(double m, double n, double k)
{
    int flat = 124, down = 47, up = 98, dive = 60, soar = 62 ;

    int desiredChar = flat*(m==n && k==m) ^
                      down*(k<m && n<k) ^
                        up*((k>m && n>k) || (n>k && k==m)) ^
                        dive*(k<m && k<n) ^
                        soar*(k>m && k>m) //^
                        //flat*(k==m || k==n)
                        ;

    return desiredChar;
}

void showHistogram(int data[], int length)
{
    int i, j, beggining, finish, height;
    beggining = smallestValue(data, length);
    finish = biggestValue(data, length);
    height = finish-beggining;
    int hist[(height+1)];
    int count[(*(&hist +1)-hist)];
    int percentages[(*(&count +1) - count)];

    //initialize histogram as every number from lowest of data to biggest of data
    for(i=0, j=beggining; i<(*(&hist +1)-hist); i++, j++)
    {
        hist[i] = j;
    }

    //initialize count as 0 for each
    for(i=0; i<(*(&count + 1)-count); i++)
    {
        count[i] = 0;
    }

    //for every element in count, assign the quantity of times an element of data appears
    for(i = 0; i<(*(&count +1)-count); i++)
    {
        for(j=0; j<length; j++)
        {
            count[i] += 1*(hist[i]==data[j]);
        }
    }
    int oneHundredPercent = length;
    //100% = length
    for(i=0; i<*(&count +1)-count;i++)
    {
        percentages[i] = count[i]*100/oneHundredPercent;
    }

    for(i=0, j=beggining; i<(*(&percentages +1)-percentages); i++, j++)
    {
        std::cout << j << '\t';
        for(int ammount=0; ammount<percentages[i]; ammount++)
        {
            std::cout << '*';
        }
        std::cout << '\n';
    }
}

void showHistogramToScale(int data[], int length)
{
    int i, j, beggining, finish, height;
    beggining = smallestValue(data, length);
    finish = biggestValue(data, length);
    height = finish-beggining;
    int hist[(height+1)];
    int count[(*(&hist +1)-hist)];
    int percentages[(*(&count +1) - count)];

    //initialize histogram as every number from lowest of data to biggest of data
    for(i=0, j=beggining; i<(*(&hist +1)-hist); i++, j++)
    {
        hist[i] = j;
    }

    //initialize count as 0 for each
    for(i=0; i<(*(&count + 1)-count); i++)
    {
        count[i] = 0;
    }

    //for every element in count, assign the quantity of times an element of data appears
    for(i = 0; i<(*(&count +1)-count); i++)
    {
        for(j=0; j<length; j++)
        {
            count[i] += 1*(hist[i]==data[j]);
        }
    }
    int oneHundredPercent = biggestValue(count, *(&count +1)-count)/2*3;
    //100% = length
    for(i=0; i<*(&count +1)-count;i++)
    {
        percentages[i] = count[i]*100/oneHundredPercent;
    }

    for(i=0, j=beggining; i<(*(&percentages +1)-percentages); i++, j++)
    {
        std::cout << j << '\t';
        for(int ammount=0; ammount<percentages[i]; ammount++)
        {
            std::cout << '+';
        }
            putchar(histchar(percentages[i-1], percentages[i+1], percentages[i]));
        std::cout << '\n';
    }
}

void showHistogramToScale(double data[], int length, double amplitude)
{
    int i, j;
    double currentInterval, smol = smallestValue(data, length), jay;

    int count[intervals(data, length, amplitude, true)];
    int percentages[(*(&count +1) - count)];


    //initialize count as 0 for each
    for(i=0; i<(*(&count + 1)-count); i++)
    {
        count[i] = 0;
    }

    //for every element in count, assign the quantity of times an element in amplitude of the interval in the dataset appears
    for(i = 0; i<length; i++)
    {
        currentInterval = (data[i] - smol) / amplitude;
        count[(int)currentInterval] += 1;
        //std::cout << '\n' << "COUNT["<< (int)currentInterval << "]: " << count[(int)currentInterval];
    }
    int oneHundredPercent = biggestValue(count, *(&count +1)-count)/2*3;
    //100% = length
    for(i=0; i<*(&count +1)-count;i++)
    {
        percentages[i] = count[i]*100/oneHundredPercent;
    }

    for(i=0, jay=smol; i<(*(&percentages +1)-percentages); i++, j++)
    {
        printf("%.3f", jay);
        std::cout << '\t' << '|';
        for(int ammount=0; ammount<percentages[i]; ammount++)
        {
            std::cout << '+';
        }
            putchar(histchar(percentages[i-1], percentages[i+1], percentages[i]));
        std::cout << '\n';
        jay += amplitude;
    }
}

int gaussianDistribution(int mu , int sigma) //sigma is the standard deviation, mu is the mean
{
    int range = sigma*2;
    int valuetoAdd = mu-range;
    int U1 = randal(0, range);
    int U2 = randal(0, range);
    //int c = U1+U2;
    //if(c<0 || c>range*2) return gaussianDistribution(mu, sigma);
    //std::cout << "VALUE TO ADD: " << valuetoAdd << '\n';
    int result = U1+U2 + valuetoAdd;

    return result;
}

void showDetails(int data[], int length)
{
    char newline = '\n';
    std::cout << "================================" << newline;
    std::cout << "mean: " << mean(data, length) << newline <<
                "median: " << median(data, length) << newline <<
                "mode: " << mode(data, length) << newline <<
                "range: " << range(data, length) << newline;
    std::cout << "-------------------------------" << newline;
    std::cout << "variance: " << variance(data, length) << newline <<
                "Standard Deviation: " << standardDeviation(data, length) << newline;
    std::cout << "================================" << std::endl;
}

void showDetails(double data[], int length)
{
    double amp = minimumAmplitude(data, length);
    char newline = '\n';
    std::cout << "================================" << newline;
    std::cout << "mean: " << mean(data, length) << newline <<
                "median: " << median(data, length) << newline <<
                "mode: " << mode(data, length, amp) << newline <<
                "range: " << range(data, length) << newline;
    std::cout << "-------------------------------" << newline;
    std::cout << "variance: " << variance(data, length) << newline <<
                "Standard Deviation: " << standardDeviation(data, length) << newline;
    std::cout << "================================" << std::endl;
}

void showDetails(double data[], int length, double amplitude)
{
    double amp = amplitude;
    char newline = '\n';
    std::cout << "================================" << newline;
    std::cout << "mean: " << mean(data, length) << newline <<
                "median: " << median(data, length) << newline <<
                "mode: " << mode(data, length, amp) << newline <<
                "range: " << range(data, length) << newline;
    std::cout << "--------------------------------" << newline;
    std::cout << "variance: " << variance(data, length) << newline <<
                "Standard Deviation: " << standardDeviation(data, length) << newline;
    std::cout << "================================" << std::endl;
}

double gaussian(double mu , double sigma) //sigma is the standard deviation, mu is the mean
{
    double U1, U2;
    U1 = randfrom(-1, 1);
    U2 = randfrom(-1, 1);
    double r = U1*U1 + U2*U2;
    if(r < (-1) || r > (1) || r==0) return gaussian(mu, sigma);
    U1 *= sigma;
    U2 *= sigma;

    double result = sqrt(-2 * log(r) / r+1*(r==0));

    return (result * U1) + mu;
}

int gaussian(int mu , int sigma) //sigma is the standard deviation, mu is the mean
{
    double U1, U2;
    U1 = randfrom(-1, 1);
    U2 = randfrom(-1, 1);
    double r = U1*U1 + U2*U2;
    if(r < (-1) || r > (1) || r==0) return gaussian(mu, sigma);
    U1 *= sigma;
    U2 *= sigma;

    double result = sqrt(-2 * log(r) / r+1*(r==0));
    return (int) ((result * U1) + mu);
}

std::pair<double, double> estimateCoefficient(int dataX[], int lengthX, int dataY[], int lengthY)
{
    std::pair<double, double> Betas;
    double meanOfX = mean(dataX, lengthX), meanOfY = mean(dataY, lengthY);
    int i;
    double residualsofX[lengthX], residualsofY[lengthY];

    for(i=0; i<lengthX && i<lengthY; i++)
    {
        residualsofX[i] = dataX[i] - meanOfX;
        residualsofY[i] = dataY[i] - meanOfY;
    }

    double sumofResiduals = 0.0;

    for(i=0; i<lengthX && i<lengthY; i++)
    {
        sumofResiduals = residualsofX[i] * residualsofY[i];
    }

    double sumoftheSquaRedresidualsofX = 0.0;

	for(i=0; i<lengthX && i<lengthY; i++)
	{
		sumoftheSquaRedresidualsofX += pow((meanOfX - dataX[i]), 2);
	}

    double Beta1 = sumofResiduals / sumoftheSquaRedresidualsofX;
    double Beta0 = meanOfY - Beta1*meanOfX;

    Betas.first = Beta0;
    Betas.second = Beta1;

    return Betas;
}

std::pair<double, double> estimateCoefficient(double dataX[], int lengthX, double dataY[], int lengthY)
{
    std::pair<double, double> Betas;
    double meanOfX = mean(dataX, lengthX), meanOfY = mean(dataY, lengthY);
    int i;
    double residualsofX[lengthX], residualsofY[lengthY];

    for(i=0; i<lengthX && i<lengthY; i++)
    {
        residualsofX[i] = dataX[i] - meanOfX;
        residualsofY[i] = dataY[i] - meanOfY;
    }

    double sumofResiduals = 0.0;

    for(i=0; i<lengthX && i<lengthY; i++)
    {
        sumofResiduals = residualsofX[i] * residualsofY[i];
    }

    double sumoftheSquaRedresidualsofX = 0.0;

	for(i=0; i<lengthX && i<lengthY; i++)
	{
		sumoftheSquaRedresidualsofX += pow((meanOfX - dataX[i]), 2);
	}

    double Beta1 = sumofResiduals / sumoftheSquaRedresidualsofX;
    double Beta0 = meanOfY - Beta1*meanOfX;

    Betas.first = Beta0;
    Betas.second = Beta1;

    return Betas;
}

std::pair<double, double> betas(double dataX[], int lengthX, double dataY[],int lengthY)
{
    double  sumOvX = 0.0, sumOfY = 0.0, crossSum = 0.0, sumovXsquared = 0.0;
    const int observations = lengthX;
    int i;
    for(i=0; i<lengthX; i++)
    {
        sumOvX += dataX[i];
        sumOfY += dataY[i];
        crossSum += dataX[i] * dataY[i];
        sumovXsquared += pow(dataX[i], 2);
        //printf("\nSumOvX: %f", sumOvX);
    }
    double Beta1 =  (observations * crossSum - sumOvX * sumOfY) /
                    (observations * sumovXsquared - pow(sumOvX, 2));
    double Beta0 =  (sumOfY * sumovXsquared - sumOvX * crossSum) /
                    (observations * sumovXsquared - pow(sumOvX, 2));
    //char n = '\n';

    //printf("Beta0:%c%.3f * %.3f - %.3f * %.3f /%c", n, sumOfY, sumovXsquared, sumOvX, crossSum, n);
    //printf("%d * %.3f - %.3f%c", observations, sumovXsquared, pow(sumOvX, 2), n);

    std::pair<double, double> Betaoneandxero;
    Betaoneandxero.first = Beta0;
    Betaoneandxero.second = Beta1;

    return Betaoneandxero;
}

std::pair<double, double> betas(int dataX[], int lengthX, int dataY[],int lengthY)
{
    double  sumOvX = 0.0, sumOfY = 0.0, crossSum = 0.0, sumovXsquared = 0.0;
    const int observations = lengthX;
    int i;
    for(i=0; i<lengthX; i++)
    {
        sumOvX += dataX[i];
        sumOfY += dataY[i];
        crossSum += dataX[i] * dataY[i];
        sumovXsquared += pow(dataX[i], 2);
        //printf("\nSumOvX: %f", sumOvX);
    }
    double Beta1 =  (observations * crossSum - sumOvX * sumOfY) /
                    (observations * sumovXsquared - pow(sumOvX, 2));
    double Beta0 =  (sumOfY * sumovXsquared - sumOvX * crossSum) /
                    (observations * sumovXsquared - pow(sumOvX, 2));
    //char n = '\n';

    //printf("Beta0:%c%.3f * %.3f - %.3f * %.3f /%c", n, sumOfY, sumovXsquared, sumOvX, crossSum, n);
    //printf("%d * %.3f - %.3f%c", observations, sumovXsquared, pow(sumOvX, 2), n);

    std::pair<double, double> Betaoneandxero;
    Betaoneandxero.first = Beta0;
    Betaoneandxero.second = Beta1;

    return Betaoneandxero;
}

double rootMeanSquaredError(int dataX[], int lengthX, int dataY[],int lengthY)
{
    std::pair<double, double> coefficients = estimateCoefficient(dataX, lengthX, dataY, lengthY);

    double predictionsOfY[lengthY];
    double Beta0 = coefficients.first, Beta1 = coefficients.second;
    int i;

    for(i=0; i<lengthX && i<lengthY; i++)
    {
        predictionsOfY[i] = Beta0 + Beta1 * dataX[i];
    }
    double meanSquared = 0.0;

    for(i=0; i<lengthX && lengthY; i++)
    {
        meanSquared += pow(predictionsOfY[i] - dataY[i], 2);
    }
    double res = meanSquared/lengthY;

    return sqrt(res);
}

double rootMeanSquaredError(double dataX[], int lengthX, double dataY[],int lengthY)
{
    std::pair<double, double> coefficients = betas(dataX, lengthX, dataY, lengthY);

    double predictionsOfY[lengthY];
    double Beta0 = coefficients.first, Beta1 = coefficients.second;
    int i;

    for(i=0; i<lengthX && i<lengthY; i++)
    {
        predictionsOfY[i] = Beta0 + Beta1 * dataX[i];
    }
    double meanSquared = 0.0;

    for(i=0; i<lengthX && lengthY; i++)
    {
        meanSquared += pow(predictionsOfY[i] - dataY[i], 2);
    }
    double res = meanSquared/lengthY;

    return sqrt(res);
}

void predictionEquation(int dataX[], int lengthX, int dataY[],int lengthY)
{
    std::pair<double, double> coefficients = betas(dataX, lengthX, dataY, lengthY);

    printf("y' = %f + %f * x%c", coefficients.first, coefficients.second, '\n');
}

void predictionEquation(double dataX[], int lengthX, double dataY[],int lengthY)
{
    std::pair<double, double> coefficients = betas(dataX, lengthX, dataY, lengthY);

    printf("y' = %f + %f * x%c", coefficients.first, coefficients.second, '\n');
}



//Linear Regression pls learn linear regression u dumb cunt
//FUCK YOU I TOLD YA I COULD LEARN IT YOU CUNT
//fucking hate you C
//oh yeah also hate that stupid language (jk love you C)

#endif // MYSTATS_H_INCLUDED
