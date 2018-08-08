/* Shemchuk Serhii K-26, Numerical methods
 * ----------
 * Project #4 Building Interpolation Polynomial
 *
 * Function: f(x) = e^(-x^2)
 *
 * */

typedef long double ld;

const ld PI = acos(-1.0);

ld func(const ld &x) {
    return (ld) exp(-x * x);
}

ld interpol_unit(const ld &a, const ld &b, const ld &k, const ld &n) {
    ld result = (b + a) / 2 + (b - a) / 2 * (ld) cos(((k + k + 1) * PI) / (n + n + 2));
    return result;
}

// Helper coefficient for calculating polynomial
ld Q(const int &n, const int &i, const ld &x, const vector<ld> &X) {
    ld result = 1;
    for (int j = 0; j < n; j++) {
        if (j != i) {
            result *= (x - X[j]) / (X[i] - X[j]);
        }
    }
    return result;
}

// Lagrange polynomial -- O(n^2)
ld L_n_3(const int &n, const ld &x, const vector<ld> &X, const vector<ld> &Y) {
    ld result = 0;
    for (int i = 0; i < n; i++) {
        result += Y[i] * Q(n, i, x, X);
    }
    return result;
}

// Lagrange polynomial -- O(n)
ld L_n_2(const int &n, const ld &x, const vector<ld> &X, const vector<ld> &Y, const vector<ld> &A) {
    ld up = 0, down = 0;

    for (int j = 0; j < n; j++) {
        if (x != X[j])
            up += (Y[j] * A[j] / (x - X[j]));
        else
            up += Y[j];
    }

    for (int j = 0; j < n; j++) {
        if (x != X[j])
            down += (A[j] / (x - X[j]));
        else
            down += 1;
    }

    return up / down;
}

inline ld step(const ld &start, const ld &h, int i) {
    return start + h * i;
}

void Lagrange(const int &n, const int &a, const int &b) {
    ld h = ((ld) b - (ld) a) / (ld) n;

    vector<ld> chebX(n), equidistX(n);
    for (int i = 0; i < n; i++) {
        // Chebyshev nodes
        chebX[i] = interpol_unit(a, b, i, n);

        // Equidistant nodes
        equidistX[i] = step(a, h, i);
    }

    vector<ld> chebY(n), equidistY(n);
    for (int i = 0; i < n; i++) {
        // Chebyshev nodes
        chebY[i] = func(chebX[i]);

        // Equidistant nodes
        equidistY[i] = func(equidistX[i]);
    }

    const int POINTS_NUMBER_TO_COMPARE = 300;
    ld step = ((ld) b - (ld) a) / (ld) POINTS_NUMBER_TO_COMPARE;

    vector<ld> points;
    for (int i = 0; i < POINTS_NUMBER_TO_COMPARE; i++) {
        ld val = (ld) a + i * step;
        points.push_back(val);
    }

    sort(points.begin(), points.end());

    const int N = (int) points.size();

    vector<ld> extraPointsForChebY(N);
    for (int i = 0; i < N; i++) {
        // O(n^2) formula
        extraPointsForChebY[i] = L_n_3(n, points[i], chebX, chebY);
    }

    vector<ld> extraPointsForEquidistY(N);
    for (int i = 0; i < N; i++) {
        // O(n^2) formula
        extraPointsForEquidistY[i] = L_n_3(n, points[i], equidistX, equidistY);
    }
}