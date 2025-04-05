# Множество Мандельброта

## О проекте

В этом проекте изучается возможность оптимизации программ с помощью векторных регистров и операций с ними. Примером для оптимизации стала отрисовка множества Мандельброта, которая требует большого количества вычислений. В программе имеется 2 режима: графический и тестовый. Тестовый режим позволяет получить количество тактов процессора, требуемое для выполнения рассчётов. Графический режим позволяет рассмотреть множество, в нём доступны перемещение с помощью стрелок на клавиатуре, а также увеличение и уменьшение масштаба с помощью клавиш 'P' и 'O' соответственно.

![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/ab69d53dd5775ec35ab5d4515d274ce37909f673/graphics.png)

## Установка

Программа имеет 3 уровня оптимизации для архитектуры X86-64 и 2 уровня для архитектуры ARM. Уровень оптимизации меняется до момента компиляции (это может поменяться в ближайшем будущем). Для выбора уровня оптимизации в файле 'source/mandelbrot.cpp' нужно изменить константу в самом верху файла. Для архитектуры X86-64 доступны значения: RENDER_VECTOR_1, RENDER_VECTOR_4 и RENDER_VECTOR_8, а для архитектуры ARM доступны только 2 из них: RENDER_VECTOR_1 и RENDER_VECTOR_4. Значения этих констант будут объяснены в слудующем разделе. Для выбора архитектуры для которой надо скомпилировать программу нужно написать в консоли: ```make x86``` и ```make arm``` для X86-64 и ARM соответственно.

## Теоретическое введение

Для оптимизации используются SIMD наборы инструкций, которые позволяют выполнять операции с несколькими числами с плавающими точками одновременно. Для X86-64 доступны оптимизации с 128 и 256 -битными регистрами. Они хранят соответственно 4 и 8 чисел с плавающей точкой одинарной точности. Для ARM доступна только оптимизация с 4 числами. В этом и заключается суть констант из предыдущего раздела:
- RENDER_VECTOR_1 - полностью отключает оптимизацию и все операции выполняются с отдельными числами.
- RENDER_VECTOR_4 - включает оптимизацию, суть которой заключается в одновременном выполнении операций с 4 точками.
- RENDER_VECTOR_8 - включает оптимизацию аналогичную предыдущей, только данные упаковываются в пачки по 8.

## Что должно получиться в теории

Теоретически уровень оптимизации RENDER_VECTOR_4 должен давать ускорение программы в ~4 раза относительно уровня RENDER_VECTOR_1, так как каждая соответствующая операция выполняется сразу для 4 чисел. Различия могут возникать из-за небольшой разницы возможностей для аналогичных операций. Для ARM архитектуры различие теоретически должно быть больше, так как в данной версии проверка на 0 четырёх чисел выполняется с помощью записи их в массив и 4 сравнений каждого по отдельности. Это достаточно долгая операция, которая выполняется на каждой итерации.
Аналогично уровень RENDER_VECTOR_8 должен давать ускорение в ~8 раз относительно уровня RENDER_VECTOR_1.

## Методика измерений

Измерение времени работы программы - достаточно непростая задача. Искомый параметр зависит от многих факторов, который практически невозможно проконтролировать. Для более точной оценки ускорения используются интринсики для соответствующих архитектур, которые позволяют получать количество тактов, выполненных процессором в момент выполнения нашей программы. Данный параметр всё равно является достаточно случайным, для наилучшей оценки количества тактов процессора, которое затрачивается на отрисовку одного экрана подоёт метод наименьших квадратов. Он также позволяет отсечь время, которое тратится на вызов функции отрисовки и на создание переменных перед циклом.

Измерения построены так: функция отрисовки вызывается с параметром, который означает количество итераций прорисовки которое она должна выполнить. Этот параметр возрастает линейно, его возрастание задано в момент запуска программы (bin/Mandelbrot --test A B C --output 'file', здесь A - первое количество итераций отрисовки, B - шаг, C - количество шагов, 'file' - имя файла в который программа запишет все экспериментальные точки) Количество тактов процессора затраченного на вызов и отрисовку несколько раз определяется с помощью интринсиков для соотвествующей архитектуры. Далее строится зависимость количества тактов затраченных на выполнения от количества итераций отрисовки. Коэффициент угла наклона этого графика теоретически должен являться средним значением для одной отрисовки.

Использование в измерениях регистра с количеством прошедших тактов процессора позволяет получить независимость от загрузки системы.

## Полученные результаты

Полученные значения приведены [здесь](/values.md).

Данные наносим на графики, с помощью метода наименьших квадратов можно оценить погрешность, так как она достаточно случайна при таких измерениях.

### Графики зависимостей количества тактов процессора затраченного на отрисовку от количества итераций

#### Без оптимизации (x86-64)

![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/ab69d53dd5775ec35ab5d4515d274ce37909f673/x86_opt1.png)

#### Оптимизация с данными упакованными по 4 элемента (x86-64)

![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/ab69d53dd5775ec35ab5d4515d274ce37909f673/x86_opt4.png)

#### Оптимизация с данными упакованными по 8 элементов (x86-64)

![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/ab69d53dd5775ec35ab5d4515d274ce37909f673/x86_opt8.png)

#### Без оптимизации (arm)

![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/ab69d53dd5775ec35ab5d4515d274ce37909f673/arm_opt1.png)

#### Оптимизация с данными упакованными по 4 элемента (arm)

![alt text](https://raw.githubusercontent.com/artemneskorodov/Mandelbrot/ab69d53dd5775ec35ab5d4515d274ce37909f673/arm_opt4.png)

### Полученные результаты количества тактов, затрачиваемых на отрисовку одного экрана:

#### Без оптимизации (x86-64)

```math
T_{x86-1} = (189 ± 1) * 10^6
```

#### Оптимизация с данными упакованными по 4 элемента (x86-64)

```math
T_{x86-4} = (47 ± 1) * 10^6
```

#### Оптимизация с данными упакованными по 8 элементов (x86-64)

```math
T_{x86-8} = (24 ± 1) * 10^6
```

#### Без оптимизации (arm)

```math
T_{arm-1} = (53 ± 1) * 10^6
```

#### Оптимизация с данными упакованными по 4 элемента (arm)

```math
T_{arm-4} = (14 ± 1) * 10^6
```


## Итог
### Полученные отношения

```math
\frac{T_{x86-1}}{T_{x86-4}} = (4.0 ± 0.1)
```

```math
\frac{T_{x86-1}}{T_{x86-8}} = (7.9 ± 0.3)
```

```math
\frac{T_{arm-1}}{T_{arm-4}} = (3.8 ± 0.3)
```
Это в пределах погрешностей измерений совпадает с примерным предсказанным результатом. Различия могут появляться из-за разных состояний системы в момент исполнения программы. Для архитектуры ARM, как и предполагалось значение этого отношения меньше 4.

