# MACHO: Multiple Alignments Column Hits Observer
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1K7gKSyI1rx_F28ohZgmQf4m9CvR5l64d?usp=sharing)<br>
AUTHORS of Multiple Alignment Column Hits Observer (MACHO)
* Daniil Nagornyi
* Vsevolod Maslenikov
* Vitalii Gagarochkin<br>
A Python tool to compare two multiple sequence alignments (MSAs) and identify matching columns.

## Описание

MACHO (Multiple Alignment Column Hits Observer) - это утилита командной строки на языке Python, предназначенная для сравнения двух множественных выравниваний последовательностей, представленных в формате FASTA. Программа находит позиции (колонки), которые являются идентичными в обоих выравниваниях, и выводит информацию об этих совпадениях в файл формата TSV, а также предоставляет сводную статистику в стандартный вывод (STDOUT).

## Использование

Программа запускается из командной строки по следующему синтаксису:

```bash
python macho.py [-h] [-hr] [-g] sequence_1 sequence_2 out
```
или

```bash
./macho.py [-h] [-hr] [-g] sequence_1 sequence_2 out
```


## Возможности

* Сравнение двух множественных выравниваний последовательностей.
* Фильтрация колонок по максимальному количеству гэпов [-g]
* Поддержка входных файлов в формате FASTA.
* Выявление колонок, которые посимвольно идентичны в обоих выравниваниях.
* Вывод найденных совпадений в файл формата TSV.
* Два режима вывода совпадений:
    * Попарные индексы совпадающих колонок.
    * Диапазоны индексов для последовательных совпадающих колонок.
* Вывод в STDOUT длин входных выравниваний и доли совпавших колонок для каждого из них.
* Интерфейс командной строки с подробной справкой.

## Установка

Для работы MACHO требуется Python 3.

*  Клонируйте репозиторий:
    ```bash
    git clone <URLРЕПОЗИТОРИЯ>
    cd <ПАПКА_С_РЕПОЗИТОРИЕМ>
    ```
