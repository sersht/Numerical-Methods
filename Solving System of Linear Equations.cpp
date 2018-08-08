/* Shemchuk Serhii K-26, Numerical methods
 * ----------
 * Project #2 Solving System of Linear Equations
 *
 * Problem: solving matrix equation using straight and approximate methods
 *
 * Using constant 3 in constructors and cycles for comfort - in this variant we have 3x3 matrix
 *
 *              ( 2   2   1 )       (  3 )
 * Variant: А = ( 2   5  -2 )   b = ( -3 )
 *              ( 1  -2  10 )       ( 14 )
 * */

const int INF = (1 << 30);

inline double Sgn(const double &value) {
    const double EPS = 1e-6;
    return (abs(value) < EPS ? 0 : value >= EPS ? 1 : -1);
}

bool is_symmetric(const vector<vector<double> > &A) {
    for (int i = 0; i < (int) A.size(); ++i) {
        for (int j = 0; j < (int) A.size(); ++j) {
            if (A[i][j] != A[j][i]) {
                return false;
            }
        }
    }

    return true;
}

void transpose_(
        const vector<vector<double> > &input,
        vector<vector<double> > &output)
{
    for (int i = 0; i < (int) input.size(); ++i) {
        for (int j = 0; j < (int) input.size(); ++j) {
            output[i][j] = input[j][i];
        }
    }
}

void cholesky_decomposition(
        const vector<vector<double> > &A,
        vector<vector<double> > &S,
        vector<vector<double> > &D,
        vector<vector<double> > &S_T)
{
    const int n = (int) A.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            S[i][j] = 0;
            D[i][j] = 0;
        }
    }

    for (int i = 0; i < n; ++i) {
        double temp = 0;
        for (int k = 0; k <= i - 1; ++k) {
            temp += (S[k][i] * S[k][i] * D[k][k]);
        }

        D[i][i] = Sgn(A[i][i] - temp);
        S[i][i] = sqrt(abs(A[i][i] - temp));

        for (int j = i + 1; j < n; ++j) {
            temp = 0;
            for (int k = 0; k <= i - 1; ++k) {
                temp += (S[k][i] * S[k][j] * D[k][k]);
            }

            S[i][j] = (A[i][j] - temp) / D[i][i] / S[i][i];
        }
    }

    transpose_(S, S_T);
}

void matrix_multiply(
        vector<vector<double> > &result,
        const vector<vector<double> > &L,
        const vector<vector<double> > &R)
{
    const int n = 3;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < (int) R[0].size(); ++j) {
            double temp = 0;
            for (int ii = 0; ii < n; ++ii) {
                temp += (L[i][ii] * R[ii][j]);
            }

            result[i][j] = temp;
        }
    }
}

void solve_back(
        const vector<vector<double> > &A,
        vector<vector<double> > &y,
        const vector<vector<double> > &b)
{
    for (int y_i = 2; y_i >= 0; --y_i) {
        double temp = 0;
        for (int i = y_i + 1; i < 3; ++i) {
            temp += (y[i][0] * A[y_i][i]);
        }

        y[y_i][0] = (b[y_i][0] - temp) / A[y_i][y_i];
    }
}

void solve_straight(
        const vector<vector<double> > &A,
        vector<vector<double> > &x,
        const vector<vector<double> > &y)
{
    for (int x_i = 0; x_i < 3; ++x_i) {
        double temp = 0;
        for (int i = 0; i < x_i; ++i) {
            temp += (x[i][0] * A[x_i][i]);
        }

        x[x_i][0] = (y[x_i][0] - temp) / A[x_i][x_i];
    }
}

void process(
        const vector<vector<double> > &S,
        const vector<vector<double> > &D,
        const vector<vector<double> > &S_T,
        const vector<vector<double> > &b,
        vector<vector<double> > &answer)
{
    // 2 step : S_T * D * y = b
    vector<vector<double> > S_TD(3, vector<double>(3));
    matrix_multiply(S_TD, S_T, D);

    vector<vector<double> > y(3, vector<double>(1));
    solve_straight(S_TD, y, b);

    // 3 step : S * x = y
    solve_back(S, answer, y);
}

void find_inversed_matrix(
        const vector<vector<double> > &A,
        const vector<vector<double> > &S,
        const vector<vector<double> > &D,
        const vector<vector<double> > &S_T)
{
    vector<vector<double> > x_1(3, vector<double>(1));
    vector<vector<double> > E_1 = {{1},
                                   {0},
                                   {0}};
    process(S_T, D, S, E_1, x_1);

    vector<vector<double> > x_2(3, vector<double>(1));
    vector<vector<double> > E_2 = {{0},
                                   {1},
                                   {0}};
    process(S_T, D, S, E_2, x_2);

    vector<vector<double> > x_3(3, vector<double>(1));
    vector<vector<double> > E_3 = {{0},
                                   {0},
                                   {1}};
    process(S_T, D, S, E_3, x_3);

    double norm1 = -INF;
    for (int i = 0; i < 3; ++i) {
        double temp = abs(A[i][0]) + abs(A[i][1]) + abs(A[i][2]);
        norm1 = max(norm1, temp);
    }

    double norm2 = -INF;
    for (int i = 0; i < 3; ++i) {
        double temp = abs(x_1[i][0]) + abs(x_2[i][0]) + abs(x_3[i][0]);
        norm2 = max(norm2, temp);
    }

    double mu = norm2 * norm1;
}

bool converge(const vector<vector<double> > &xk, const vector<vector<double> > &xkp) {
    const double eps = 1e-6;

    double norm = 0;
    for (int i = 0; i < 3; i++) {
        norm += (xk[i][0] - xkp[i][0]) * (xk[i][0] - xkp[i][0]);
    }

    return sqrt(norm) < eps;
}

void gauss_seidel_method(
        const vector<vector<double> > &A,
        vector<vector<double> > &xx,
        const vector<vector<double> > &b)
{
    int diag_sum = 0;
    for (int i = 0; i < 3; i++) {
        diag_sum += abs(A[i][i]);
    }

    int sum = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }
    }

    if (diag_sum < sum) {
        cout << "Не гарантується, що метод Зейделя зійдеться (не виконуються достатні умови)" << endl;
    }

    vector<vector<double> > p(3, vector<double>(1));
    int count = 0;
    do {
        for (int i = 0; i < 3; i++)
            p[i][0] = xx[i][0];

        for (int i = 0; i < 3; i++) {
            double var = 0;

            for (int j = 0; j < i; j++)
                var += (A[i][j] * xx[j][0]);

            for (int j = i + 1; j < 3; j++)
                var += (A[i][j] * p[j][0]);

            xx[i][0] = (b[i][0] - var) / A[i][i];
        }

        count++;
    } while (!converge(xx, p));
}

void solve_(
        const vector<vector<double> > &A,
        const vector<vector<double> > &b)
{
    vector<vector<double> >
            S(3, vector<double>(3)),
            D(3, vector<double>(3)),
            S_T(3, vector<double>(3));

    if (!is_symmetric(A)) {
        cout << "Несиметрична!\n";
        return;
    }

    // 1 step : A = S_T * D * S
    cholesky_decomposition(A, S, D, S_T);

    vector<vector<double> > x(3, vector<double>(1));
    process(S, D, S_T, b, x);

    // Answer output
    cout << "Розв'язок рівняння [A](3, 3) * [x](3, 1) = [b](3, 1) прямим методом :\n";
    for (int i = 0; i < 3; ++i) {
        cout << x[i][0] << "\n";
    }

    // Check if the answer is correct
    cout << "\n" << "Перевірка (помножимо [A](3, 3) на отриманий розв'язок) :\n";
    vector<vector<double> > res(3, vector<double>(1));
    matrix_multiply(res, A, x);
    for (int i = 0; i < 3; ++i) {
        cout << res[i][0] << "\n";
    }

    // Calculating the determinant of the matrix
    double det = 1;
    for (int i = 0; i < 3; ++i) {
        det *= (D[i][i] * S[i][i] * S[i][i]);
    }
    cout << "\n" << "Det(A) = " << det << "\n";

    // Calculating A^(-1)
    find_inversed_matrix(A, S_T, D, S);

    // Calculating answer
    vector<vector<double> > _x_(3, vector<double>(1, 0));
    gauss_seidel_method(A, _x_, b);
}