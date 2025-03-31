# Множество Мандельброта

## О проекте

В этом проекте проведено изучение возможностей для оптимизации кода на языках C/C++ с использованием упакованных данных в xmm/ymm регистрах. Множество Мандельброта используется как пример для сложных и долгих вычислений, которые требуют много итераций циклов.

## В проекте использовались

- Библиотека SFML для графического отображения
- Функция clock_gettime() с параметром CLOCK_PROCESS_CPUTIME_ID, который позволяет получить время работы нашего процесса
- Сайт intel intrinsics guide, который позволил находить нужные интринсики для использования

## Что должно получиться в теории

Для оптимизации использовались инструкции, позволяющие делать вычисления не для одной точки, а для 4 или 8. Теоретически это должно ускорять работу программы в ~4 и ~8 раз соотственно. Для измерения времени использовался метод наименьших квадратов. Это позволяет отбросить время которое тратится на задание констант и вызов функций. Производятся измерения зависимости времени исполнения от количества итераций, соответственно коэффициент наклона полученной прямой должен являться временем затрачиваемым на рассчёты для рендера одного экрана.

## Полученные результаты

Полученные значения приведены здесь.

Данные наносим на графики, с помощью метода наименьших квадратов можно оценить погрешность, так как её можно считать достаточно случайной.

### Без оптимизации
График зависимости времени рендеринга от количества итераций при отсутствии оптимизации:
![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/3dd94e5f8ad5876d68269f2d1b6581e864cfda2e/img/desmos-graph.png.png)
Время затраченное на просчёт одного кадра: T1 = (0.074956 ± 0.000099) секунд
### Оптимизация с данными упакованными по 4 элемента
График зависимости времени рендеринга от количества итериций при оптимизации с использованием 4-ёх упакованных элементов
![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/3dd94e5f8ad5876d68269f2d1b6581e864cfda2e/img/desmos-graph-2.png.png)
Время затраченное на просчёт одного кадра: T2 = (0.018268 ± 0.000020) секунд
### Оптимизация с данными упакованными по 8 элементов
График зависимости времени рендеринга от количества итериций при оптимизации с использованием 8-ми упакованных элементов
![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/3dd94e5f8ad5876d68269f2d1b6581e864cfda2e/img/desmos-graph-3.png.png)
Время затраченное на просчёт одного кадра: T3 = (0.009500 ± 0.000011) секунд


## Итог
Полученные отношения T1 / T2 = 4,10; T1 / T3 = 7,89. Это в пределах погрешностей измерений совпадает с примерным предсказанным результатом. Различия могут появляться из-за разных состояний системы в момент исполнения программы.


