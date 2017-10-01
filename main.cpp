#include <cstdio>
#include <fstream>
#include <cmath>
#include <streambuf>
#include <string>
#include <list>

using namespace std;

double** parsePoints(list<string> strPoints)
{
    int dim1Length = strPoints.size(); //Antalet noder
    int dim2Length = 2;
    double** points = new double*[dim1Length];
    int nodeIndex = 0;
    int pointIndex = 0;
    string strNodes[2];
    for (string strPoint : strPoints)
    {
        strNodes[0].clear(); //Tried clearing entire array but didn't work
        strNodes[1].clear();

        for (char& c : strPoint)
        {
            if (c == ',')
            {
                nodeIndex = 1;
            }
            else
            {
                strNodes[nodeIndex] += c;
            }
        }
        double* node = new double[dim2Length];
        node[0] = stod(strNodes[0]);
        node[1] = stod(strNodes[1]);

        points[pointIndex++] = node;

        nodeIndex = 0;
    }

    return points;
}

double** csvToMultiDim(string fileName, int *nrP)
{
    ifstream t(fileName);
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

    *nrP = unprocessedPoints.size();

    return parsePoints(unprocessedPoints);
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
    double *valueList = new double[2];
    valueList[0] = newB;
    valueList[1] = newM;

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

        delete[] sg;
    }

    double* valueList = new double[2];
    valueList[0] = b;
    valueList[1] = m;

    return valueList;
}

int main()
{
    string csvPath = "./data.csv";
    int nrPoints;
    double **points = csvToMultiDim(csvPath, &nrPoints);
    double learningRate = 0.0001;
    double initialB = 0;
    double initialM = 0;
    int numIterations = 1000;

    printf("Starting gradient descent at b = %f, m = %f, error = %f\n", initialB, initialM, computeErrorForLineGivenPoints(initialB, initialM, points, nrPoints));
    printf("Running...\n");

    double *finalVars = gradientDescentRunner(points, nrPoints, initialB, initialM, learningRate, numIterations);
    
    printf("After %d iterations on %d nodes: b = %f, m = %f, error = %f\n", numIterations, nrPoints, finalVars[0], finalVars[1], computeErrorForLineGivenPoints(finalVars[0], finalVars[1], points, nrPoints));
    
    for (int i = 0; i < nrPoints; i++)
        delete[] points[i];

    delete[] points;

    return 0;
}