/* Shemchuk Serhii K-26, Numerical methods
 * ----------
 * Project #3 Solving System of Nonlinear Equations
 *
 * Fixed z - third variable
 *
 * */

typedef long double ld;

ld Xbegin(const ld &Z) {
    return (ld) sqrt(((ld) 328 - ((ld) 9 * Z + (ld) 2) * ((ld) 9 * Z + (ld) 2)) / (ld) 315);
}

ld YBegin(const ld &X, const ld &Z) {
    return phi2(X, Z);
}

ld phi1(const ld &Y, const ld &Z) {
    return (ld) 2 * (ld) sqrt(Z - Y * Y);
}

ld phi2(const ld &X, const ld &Z) {
    return (ld) 3 * (ld) sqrt((ld) 1 - X * X - Z * Z / (ld) 4);
}

ld dPhi1_Dy(const ld &Y, const ld &Z) {
    return (ld) 2 * Y / (ld) sqrt(Z - Y * Y);
}

ld dPhi2_Dx(const ld &X, const ld &Z) {
    return (ld) 3 / (ld) sqrt((ld) 1 - X * X - Z * Z / (ld) 4);
}

ld PHI(ld &x, ld &y, const ld &z, const ld &T) {
    ld newx = x - T * (x * x / (ld) 4 + y * y - z);
    x = newx;

    ld newy = y - T * (x * x + y * y / (ld) 9 + z * z / (ld) 4 - 1);
    y = newy;
}

ld Tau(const ld &X, const ld &Y, const ld &Z, const ld &AREA) {
    ld val1 = max((ld) 3 / (ld) 2, max((ld) 3 / (ld) 2 * (X - AREA), (ld) 3 / (ld) 2 * (X + AREA)));
    ld val2 = max((ld) 20 / (ld) 9, max((ld) 20 / (ld) 9 * phi2(X - AREA, Z), (ld) 20 / (ld) 9 * phi2(X + AREA, Z)));
    return (ld) 2 / max(val1, val2);
}

bool check(const ld &Z) {
    return
        !((ld) 2 * ((ld) sqrt((ld) 17) - (ld) 4) < Z &&
        Z < (ld) 2 / (ld) 9 * ((ld) sqrt((ld) 82) - (ld) 1));
}

void Run() {
    const ld AREA = 0.01;

    // Fixing Z-variable
    cout << "Введіть змінну Z (зафіксуйте третій вимір):\n";

    ld Z;
    while (1) {
        cin >> Z;
        if (check(Z))
            cout << "При заданому Z рішень даної системи не існує. Введіть ще раз.\n";
        else
            break;
    }

    // check Z in (2 * (sqrt(17) - 4) ; 2 / 9 * (sqrt(82) - 1) )
    // обработать случай когда один из эллипсов - мнимый, когда нет точек пересечения

    // Set precision
    cout << "\nВведіть точність наближення методу:\n";
    ld EPS;
    cin >> EPS;

    cout << "\nВведіть точність обрання початкового наближення (від цього залежить збіжність)\n"
            " як кількість знаків після коми (не більше 6):\n";

    int st;
    cin >> st;
    while (st < 1 || st > 6) {
        cout << "\nВведіть точність обрання початкового наближення (від цього залежить збіжність)\n"
                " як кількість знаків після коми (не більше 6):\n";
        cin >> st;
    }

    // Checking the convergence conditions
    ld x = Xbegin(Z);
    x = (ld) round(x * pow(10, st)) / pow(10, st);

    ld y = YBegin(x, Z);
    y = (ld) round(y * pow(10, st)) / pow(10, st);

    ld T = Tau(x, y, Z, AREA);

    ld x_left = x - AREA, x_right = x + AREA;
    ld y_left = y - AREA, y_right = y + AREA;

    ld checking = norm(PHI(x_left, x_right, Z, T) - PHI(y_left, y_right, Z, T)) / norm(x - y);
    if ((ld) 1e-6 < checking && checking + (ld) 1e-6 < 1) {
        cout << "Не виконується умова збіжності!";
        return;
    }

    cout << "\nОбласть початкового наближення обрана як: \n";
    cout << "x : |x - " << x << "| <= " << AREA << endl;
    cout << "y : |y - " << y << "| <= " << AREA << endl << endl;


    // Calculating answer
    ld old_x = x, old_y = y;

    int iter_cnt = 0;
    do {
        iter_cnt++;

        old_x = x;
        old_y = y;

        PHI(x, y, Z, T);

    } while (abs(x - old_x) < EPS && abs(y - old_y) < EPS);
}