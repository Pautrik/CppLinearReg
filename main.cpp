#include <cstdio>
#include <fstream>
#include <cmath>
#include <streambuf>
#include <string>

using namespace std;

double** csvToArr(string fileName, int *nrP)
{
    ifstream t(filename);
    string str((istreambuf_iterator<char>(t)), (istreambuf_iterator<char>()));
    list<string> unprocessedPoints;
    
    string tempStr = "";
    for(char& c : str) 
    {
        if (c == '\n')
        {
            unprocessedPoints.push_back(tempStr);
            tempStr = "";
        }
        else
            tempStr += c;
    }
}

double** parsePoints(list<string> strPoints)
{
    double** points = new double*[];
    string strNodes[] = {"", ""};
    int nodeIndex = 0;
    for (string strPoint : strPoints)
    {
        for (char& c : strPoint)
        {
            if (c == ',')
            nodeIndex = 1;
            else
            strNodes[nodeIndex] += c;
        }
        nodeIndex = 0;
    }
}

double computeErrorForLineGivenPoints(double b, double m, double **points, int nrPoints)
{
    int totalError = 0;
    for (int i = 0; i < nrPoints; i++)
    {
        double x = points[i][0];
        double y = points[i][1];
        totalError += pow(y - (m * x + b), 2);
    }
    return totalError;
}

double* stepGradient(double bCurrent, double mCurrent, double **points, int nrPoints, double learningRate) {
    double bGradient = 0;
    double mGradient = 0;
    double N = nrPoints;
    
    for (int i = 0; i < nrPoints; i++)
    {
        double x = points[i][0];
        double y = points[i][1];
        bGradient += -(2/N) * (y - ((mCurrent * x) + bCurrent));
        mGradient += -(2/N) * x * (y - ((mCurrent * x) + bCurrent));
    }
    double newB = bCurrent - learningRate * bGradient;
    double newM = mCurrent - learningRate * mGradient;
    double valueList[2] = { newB, newM };

    return valueList;
}

double* gradientDescentRunner(double **points, int nrPoints, double startingB, double startingM, double learningRate, int numIterations) {
    double b = startingB;
    double m = startingM;
    for (int i = 0; i < numIterations; i++)
    {
        double *sg = stepGradient(b, m, points, nrPoints, learningRate);
        b = sg[0];
        m = sg[1];
    }
    double valueList[2] = { b, m };

    return valueList;
}

int main()
{
    string csvPath = "./data.csv";
    int *nrP;
    double **points = csvToArr(csvPath, nrP);
    int nrPoints = *nrP;
    double learningRate = 0.0001;
    double initialB = 0;
    double initialM = 0;
    int numIterations = 1000;

    printf("Starting gradient descent at b = %f, m = %f, error = %f", initialB, initialM, computeErrorForLineGivenPoints(initialB, initialM, points, nrPoints));
    printf("Running...");

    double *finalVars = gradientDescentRunner(points, nrPoints, initialB, initialM, learningRate, numIterations);
    
    printf("After %d iterations b = %f, m = %f, error = %f", numIterations, finalVars[0], finalVars[1], computeErrorForLineGivenPoints(finalVars[0], finalVars[1], points, nrPoints));

    return 0;
}