#include "nr3.h"
#include "banded.h"

int main()
{
    int i, j;

    VecDoub x(7), b(7);
    MatDoub a(7, 4);
    

    a[0][0] = 0.0000000; // arbitrary: not used
    a[0][1] = 0.0000000; // arbitrary: not used
    a[0][2] = 0.5297000;
    a[0][3] = 0.9304360;
    a[1][0] = 0.0000000; // arbitrary: not used
    a[1][1] = 0.0668422;
    a[1][2] = 0.7226600;
    a[1][3] = 0.6711490;
    a[2][0] = 0.6316350;
    a[2][1] = 0.8847070;
    a[2][2] = 0.5194160;
    a[2][3] = 0.6515190;
    a[3][0] = 0.2624530;
    a[3][1] = 0.7621980;
    a[3][2] = 0.7533560;
    a[3][3] = 0.9092080;
    a[4][0] = 0.2727100;
    a[4][1] = 0.8976560;
    a[4][2] = 0.2749070;
    a[4][3] = 0.5162920;
    a[5][0] = 0.2470390;
    a[5][1] = 0.4865170;
    a[5][2] = 0.8461670;
    a[5][3] = 0.8309650;
    a[6][0] = 0.9910370;
    a[6][1] = 0.6792960;
    a[6][2] = 0.7664950;
    a[6][3] = 0.0000000; // arbitrary: not used

    x[0]    = 0.4159990;
    x[1]    = 0.3835020;
    x[2]    = 0.3834160;
    x[3]    = 0.2377740;
    x[4]    = 0.0726859;
    x[5]    = 0.3592650;
    x[6]    = 0.0345721;

    cout << "Banded matrix [a]:" << endl;
    cout << fixed << setprecision(7);
    for (i = 0; i < a.nrows(); i++) {
        for (j = 0; j < a.ncols(); j++) {
            cout << setw(12) << a[i][j];
        }
        cout << endl;
    }
    cout << endl;

    cout << "Vector [x]:" << endl;
    for (i = 0; i < x.size(); i++) {
        cout << setw(18) << x[i] << endl;
    }
    cout << endl;

    // multiply a times x, giving b
    banmul(a, 2, 1, x, b);

    cout << "After multiplying [a] times [x] this is [b]:" << endl;

    for (i = 0; i < b.size(); i++) {
      cout << setw(18) << b[i] << endl;
    }
    cout << endl;

    //
    // First, save original value of x
    //
    VecDoub xsave(x);
    //
    // The constructor does the decomposition
    //
    Bandec banded(a, 2, 1);

    //
    // Now solve for x
    //
    banded.solve(b, x);

    cout << "After solving [a] times [x] = [b] with bandec.solve():" << endl;
    cout << "    Original     Solved" << endl 
         << "       x           x" << endl;
    for (i = 0; i < x.size(); i++) {
      cout << setw(12) << xsave[i] << setw(12) << x[i] << endl;
    }

    return 0;
}
