# Metagenomic analysis: data processing, taxonomic profiling and evaluation
## Вступление:
В данном задании я буду проводить анализ секвенирования микробного сообщества, который включает в себя следующие шаги: 
1. Проведение оценки качества данных и их очистка.
2. Проведение обработки данных с использованием пайплана [dada2](https://benjjneb.github.io/dada2/tutorial_1_8.html) для получения таблицы ASV на уровне родов и видов.
3. Анализ полученных таблиц:
   - рассчет индексов биоразнообразия, сравнение результатов между группами образцов
   - анализ на специфические таксоны (LefSE) с использованием платформы [microbiom](https://www.microbiomeanalyst.ca/)
4. Формирование гипотезы о типе отсеквенированного микробиомного сообщества

Для анализа я выбрала набор данных №3.

## Предобработка данных
Прежде всего необходимо скачать данные с гугл диска и создать необходимые папки:
```{bash}
conda activate bi_env
mkdir metagenomics_SPBU
cd metagenomics_SPBU/
mkdir Task3  # сюда перемещаем скачанные последовательности
cd Task3
ls 
```
<img width="1000" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/7358be88-f377-48e6-ac8d-b8397310e7da">

Посмотрим на качество наших данных с помощью [fasQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```{bash}
mkdir fastqc_res1
fastqc -o fastqc_res1 *fastq.gz
cd fastqc_res1/
```

После визуального анализа полученных графиков можно сделать вывод:
1) качество прямых ридов падает к концу до значения чуть менее 20
2) качество обратных ридов страдает в некоторых случаях в начале рида, а также стабильно ухудшается (при чем для большей части ридов на подходе к концу)
3) Per base sequence content не соответствует обычному виду: очень большие скачки в графике для соседних позиций оснований
4) GC-контент представлен распределенияем с несколькими пиками, что может свидетельствовать о каких-то контаменциях образцов
5) в обратных ридах часто встречаются на конце неопределенные основания
6) не обнаружено адаптеров
