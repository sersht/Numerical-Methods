/* Shemchuk Serhii K-26, Numerical methods
 * ----------
 * Project #1 Finding roots of Nonlinear Equation
 *
 * Problem: searching for the roots of the nonlinear equation using approximate methods
 *
 * Specification of methods: relaxation and modified Newton
 *
 * Function: f(x) = x^3 + 4 * sin(x)
 *
 * */

// Function
inline double f(const double &x) {
    return (x * x * x + 4 * sin(x));
}

// Derivative
inline double df(const double &x) {
    return (3 * x * x + 4 * cos(x));
}

// Second derivative
inline double ddf(const double &x) {
    return (6 * x - 4 * sin(x));
}

void Relaxation(const double &left, const double &right) {
    // Output precision
    cout << "Введіть кількість знаків після коми які ви хочете побачити\n";

    int n = 3;
    cin >> n;
    n = max(0, n);

    // Entering precision
    double EPS = 1e-9;
    cout << "Введіть '+' якщо бажаєте точність задати самостійно, або '-' якщо залишити за замовчуванням\n";

    char ans;
    cin >> ans;
    if (ans == '+') {
        cout << "Введіть точність Eps\n";
        cin >> EPS;
    }
    while (EPS >= 1) {
        cout << "Замала точність - введіть ще раз, або вийдіть з програми\n";
        cin >> EPS;
    }

    // Entering start point
    double start = right;
    cout << "Введіть '+', якщо бажаєте початкове наближення та '-' якщо залишити за замовчуванням\n";
    cin >> ans;
    while (ans != '+' && ans != '-') {
        cout << "Неправильний формат вхідних даних, спробуйте ще\n";
        cin >> ans;
    }
    if (ans == '+') {
        cin >> start;
        while (start < left || start > right) {
            cout << "Вийшли за межі границь\n";
            cout << start << "\n";
            cin >> start;
        }
    }

    cout << "Початкове наближення = " << fixed << setprecision(n) << start << "\n";

    // The conditions for convergence were checked in the report
    const double m = (4);
    const double M = (3 + cos(1.0));

    const double c = 2 / (m + M); // costant needed for this method
    double prev = 0; // x[n - 1]
    double next = start; // x[n]

    // Lower bound of iterations number
    cout << "\nНижня границя оцінюваної кількості ітерацій : ";

    const double q = (M - m) / (M + m);
    const double low = log(1.0 / abs(q));
    const double high = log(abs(start) / EPS);

    int iter_count_lower_bound = int(high / low) + 1;
    cout << iter_count_lower_bound << "\n\n";

    // Calculating answer
    int iter_count = 0;
    while (abs(next - prev) > EPS) {
        ++iter_count;

        prev = next;
        next = prev - c * f(prev); // main step

        cout << "Ітерація № " << iter_count << " = " << fixed << setprecision(n) << next << "\n";
    }

    // Finded answer
    cout << "Знайдений корінь = " << fixed << setprecision(n) << next << '\n';
}

void ModifiedNewton() {
    const int MAX_POSSIBLE = 10000;
    const double INF = 1e6;

    // Output precision
    int n = 12;
    cout << "Введіть кількість знаків після коми які ви хочете побачити\n";
    cin >> n;
    n = max(0, n);

    // Entering precision
    double EPS = 1e-9;
    cout << "Введіть '+' якщо бажаєте точність задати самостійно, або '-' якщо залишити за замовчуванням\n";

    char ans;
    cin >> ans;
    if (ans == '+') {
        cout << "Введіть точність Eps\n";
        cin >> EPS;
    }

    // Entering start point
    cout << "Введіть границі (послідовно ліву та праву): \n";

    int l, r;
    cin >> l >> r;
    double start = (r + l) / 2;

    // Checking the convergence conditions
    // Функція на відрізку перетинається з нулем, друга похідна не змінює знак на цьому відрізку
    // та будь-якого x, що є [l, r] та викон. умова f(x) * ddf(x) > 0
    if (f(l) * f(r) >= 0) {
        cout << "Метод не зійдеться!\n";
        return;
    }

    if (f(start) * ddf(start) <= 0) {
        cout << "Не гарантується що метод зійдеться\n";
    }

    cout << "Початкове наближення = " << fixed << setprecision(n) << start << "\n";

    // Lower bound of iterations number
    cout << "\nНижня границя оцінюваної кількості ітерацій : ";
    const double m1 = 4; // min(df) on [-1; 1]
    const double M2 = ddf(1); // max(ddf) on [-1; 1]

    const double q = (M2 * abs(start)) / (2 * m1);
    const double low = log(1.0 / abs(q));
    const double high = log(abs(start) / EPS);

    double iter_count_lower_bound = log((abs(high) / abs(low)) + 1) + 1;
    cout << iter_count_lower_bound << "\n\n";

    // Calculating answer
    const double DF = df(start);
    double prev;
    double next = start;

    int iter_count = 0;
    while (abs(f(next)) > EPS) {
        iter_count++;

        prev = next;
        next = prev - f(prev) / DF;

        cout << "Ітерація № " << iter_count << " = " << fixed << setprecision(n) << next << "\n";

        if (iter_count >= MAX_POSSIBLE) {
            cout << "\nТут щось не те, забагато ітерацій, я не впевнена що метод зараз зійдеться...\n";
            return;
        }
    }

    // Finded answer
    cout << "Знайдений корінь = " << fixed << setprecision(n) << next << '\n';
}