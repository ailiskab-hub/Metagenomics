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

Полученные файлы сохранены в папке [fastqc_res1](https://github.com/ailiskab-hub/Metagenomics/tree/main/fastqc_res1)

После визуального анализа полученных графиков можно сделать вывод:
1) качество прямых ридов падает к концу до значения чуть менее 20
2) качество обратных ридов страдает в некоторых случаях в начале рида, а также стабильно ухудшается (при чем для большей части ридов на подходе к концу)
3) Per base sequence content не соответствует обычному виду: очень большие скачки в графике для соседних позиций оснований
4) GC-контент представлен распределенияем с несколькими пиками, что может свидетельствовать о каких-то контаменциях образцов
5) в обратных ридах часто встречаются на конце неопределенные основания
6) не обнаружено адаптеров

Далее необходимо провести обрезку ридов, чтобы не работать с некачественными данными. Я буду использовать [TrimmomaticPE](http://www.usadellab.org/cms/?page=trimmomatic) для обрески спаренных ридов
Зададим следующие параметры обрезки:
- LEADING:22
- TRAILING:22
- SLIDINGWINDOW:4:24
- MINLEN:36
Для удобства я буду использовать snakemake. [Snakefile](https://github.com/ailiskab-hub/Metagenomics/blob/main/Snakefile) с командой находится в репозитории

Создадим папку `trimmed` для хранения обрезанных последовательностей

Далее проведем обрезку последовательностей:
```{bash}
cp /home/alisa/metagenomics_SPBU/Task3/*gz .
snakemake --cores=all -p
rm *_001* #убираем необрезанные файлы из папки
mkdir unpair_trim/
mv *unpaired* unpair_trim/ #перемещаем неспаренные риды в отдельную папочку
```
Посмотрим на результаты тримминга:
```{bash}
Input Read Pairs: 45330 Both Surviving: 39229 (86.54%) Forward Only Surviving: 3642 (8.03%) Reverse Only Surviving: 955 (2.11%) Dropped: 1504 (3.32%)

Input Read Pairs: 42061 Both Surviving: 38732 (92.09%) Forward Only Surviving: 1428 (3.40%) Reverse Only Surviving: 1138 (2.71%) Dropped: 763 (1.81%)

Input Read Pairs: 50056 Both Surviving: 42448 (84.80%) Forward Only Surviving: 5430 (10.85%) Reverse Only Surviving: 585 (1.17%) Dropped: 1593 (3.18%)

Input Read Pairs: 52523 Both Surviving: 44773 (85.24%) Forward Only Surviving: 5155 (9.81%) Reverse Only Surviving: 838 (1.60%) Dropped: 1757 (3.35%)

Input Read Pairs: 57063 Both Surviving: 51964 (91.06%) Forward Only Surviving: 2309 (4.05%) Reverse Only Surviving: 1529 (2.68%) Dropped: 1261 (2.21%)

Input Read Pairs: 54928 Both Surviving: 46630 (84.89%) Forward Only Surviving: 5796 (10.55%) Reverse Only Surviving: 680 (1.24%) Dropped: 1822 (3.32%)

Input Read Pairs: 56311 Both Surviving: 48972 (86.97%) Forward Only Surviving: 4632 (8.23%) Reverse Only Surviving: 1004 (1.78%) Dropped: 1703 (3.02%)

Input Read Pairs: 53322 Both Surviving: 47808 (89.66%) Forward Only Surviving: 3056 (5.73%) Reverse Only Surviving: 1138 (2.13%) Dropped: 1320 (2.48%)

Input Read Pairs: 61817 Both Surviving: 55815 (90.29%) Forward Only Surviving: 3182 (5.15%) Reverse Only Surviving: 1256 (2.03%) Dropped: 1564 (2.53%)

Input Read Pairs: 67072 Both Surviving: 59660 (88.95%) Forward Only Surviving: 3875 (5.78%) Reverse Only Surviving: 1785 (2.66%) Dropped: 1752 (2.61%)

Input Read Pairs: 74203 Both Surviving: 68186 (91.89%) Forward Only Surviving: 2927 (3.94%) Reverse Only Surviving: 1556 (2.10%) Dropped: 1534 (2.07%)

Input Read Pairs: 74417 Both Surviving: 68317 (91.80%) Forward Only Surviving: 3127 (4.20%) Reverse Only Surviving: 1584 (2.13%) Dropped: 1389 (1.87%)
```

В целом, сохранен достаточно большой процент ридов (более 90% почти в каждом образце)

Далее снова посмотрим на качество последовательностей посде тримминга. Будем смотреть именно на спаренные последовательности
```{bash}
cd ..
mkdir fastqc_res2
fastqc -o fastqc_res2 ./trimmed/*fastq.gz
```
Полученные файлы сохранены в папке [fastqc_res2](https://github.com/ailiskab-hub/Metagenomics/tree/main/fastqc_res2)

После визуального анализа полученных графиков можно сделать вывод:
1) качество прямых и обратных ридов немного падает к концу, но в целом соответствует значению более 20
2) Per base sequence content не соответствует обычному виду: очень большие скачки в графике для соседних позиций оснований
3) GC-контент представлен распределенияем с несколькими пиками
4) в обратных ридах на концах встречаются неопределенные основания

Распакуем обрезанные последовательности `gunzip trimmed/*gz`

## Обработка данных

Далее работа будет проходить по пайплайну [dada2](https://benjjneb.github.io/dada2/tutorial_1_8.html), в [скрипте на R](https://github.com/ailiskab-hub/Metagenomics/tree/main/R_script)

В результате мы получили ASV table и результаты таксономического анализа. Также я создала файл `metadata.csv`, в котором хранятся данные образцов. Все это ледит в папке [`analysis`](https://github.com/ailiskab-hub/Metagenomics/tree/main/analysis)


## Анализ полученных таблиц
Чтобы проанализировать полученные результаты я отправляюсь в раздел `Marker Data Profiling` на сайте [microbiom](https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/upload/OtuUploadView.xhtml)

### Community Profiling
**Rarefaction Curve Analysis**

Rarefaction Curve Analysis - помогает понять, насколько полно исследовано сообщество и как изменяется количество обнаруженных видов при увеличении числа семплов.

По оси X отображается количество семплов, а по оси Y — количество уникальных видов. График позволяет оценить, насколько сообщество насыщено видами.

<img width="673" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/e4965b54-de41-457d-86d5-114aac606feb">

Анализ кривой: Кривая стабилизируется и выходит на плато достаточно быстро, это может указывать на то, что для полного охвата разнообразия в данном сообществе достаточно имеющихся данных.

**Alpha diversity**

<img width="657" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/2b65abb7-19f5-42a1-a636-0fd3043c5db2">

<img width="637" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/c102b964-da8b-4100-99b4-250b5ad46f04">

<img width="642" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/1a5d63ee-c19d-4365-8574-63b4614d236d">

Результаты, посчитанные мной в R:

```{R}
   Wilcoxon rank sum exact test

data:  rich$Chao1[1:6] and rich$Chao1[7:12]
W = 15, p-value = 0.6991
alternative hypothesis: true location shift is not equal to 0


	Wilcoxon rank sum exact test

data:  rich$Shannon[1:6] and rich$Shannon[7:12]
W = 14, p-value = 0.5887
alternative hypothesis: true location shift is not equal to 0


	Wilcoxon rank sum test with continuity correction

data:  rich$Simpson[1:6] and rich$Simpson[7:12]
W = 15.5, p-value = 0.7457
alternative hypothesis: true location shift is not equal to 0
```
Основываясь на полученных данных можно утверждать что различия в альфа-разнообразии между контрольными и тестовыми образцами не несут статистической значимости.

_Для расчета статистики я использовала U-критерий Манна-Уитни,так как в выборках всего по 6 наблюдений, поэтому лучше использовать непараметрический критерий_

**Beta diversity**

<img width="593" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/904ada86-d18d-4def-9577-6d45636ea434">

<img width="601" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/006aa006-f11a-410c-aa64-5f6289a3f57f">


График по семействам:

<img width="480" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/82709ace-0a2f-4cb4-8dad-307c0e938224">

Основываясь на полученных данных, можно сказать, что различия в обилии микроорганизмов между контролем и тестом не выявлено. 

**Core Microbiome**<a name="core"></a>

Данный анализ позволяет определтить ключевые микроорганизмы, стабильно присутствующие в образцах. Это понадобится для дальнейшего определения микробного сообщества.

<img width="374" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/9f1d224d-442b-46a7-9355-28b8784c2304">


### Comparison & Classification
**Linear Discriminant Analysis Effect Size (LEfSe)**

Linear Discriminant Analysis Effect Size (LEfSe) — это метод для выявления статистически значимых различий между группами. Используется для выявления таксонов, которые значительно различаются между разными группами.

<img width="561" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/dc30de33-ceb2-4a71-baee-0629b1b4b42d">

<img width="397" alt="image" src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/bc44cd99-32f6-4a5e-860b-e61823efae26">

Таким образом, в наших контрольных данных присутствует род Bifidobacterium и Ruminococcus, а в тестовых - Catenibacterium и Granulicatella. 
И эти различия являются статистически значимыми.

## Формирование гипоотезы
Чтобы выяснить какое микробное сообщество было отсеквенировано нужно обратиться к данным, полученным из анализа [`Core Microbiom`](#core)

Данные в виде [таблицы](https://github.com/ailiskab-hub/Metagenomics/blob/main/core_microbiome.csv) представлены в репозитории.

1) _Faecalibacterium sp._
   
   У людей по всему миру широко распространены _Faecalibacterium_ - они были обнаружены в 85% образцах микробиоты кишечника, а представители рода _Faecalibacterium_ считаются повсеместно распространенными в желудочно-кишечном тракте здоровых людей [1]

2) _Collinsella and Bifidobacterium_

   _Collinsella and Bifidobacterium_ являются важными бактериями в развитии микробиоты кишечника. У пожилых людей по мере старения наблюдается снижение уровня популяций бифидобактерий в кишечнике [2]

3) _Streptococcus_

   Род _Streptococcus_ включает около 100 видов грамположительных бактерий, обитающих в различных средах и имеющих множество различий. Одни стрептококки патогенны, а другие являются нормальными обитателями полости рта и желудочно-кишечного тракта [3]

Судя по перечисленным родам, сообщество представляет собой микробиом кишечника. Такие роды, как _Faecalibacterium, Bifidobacterium, Collinsella, Streptococcus, Prevotella_ и другие, обычно ассоциируются с микробиотой кишечника. Присутствие этих родов позволяет предположить, что данное сообщество может быть связано с желудочно-кишечным трактом.

<img
width="641"
alt="image"
src="https://github.com/ailiskab-hub/Metagenomics/assets/112699940/461ff9d3-6fdc-44aa-9bf9-c8d06857aec2">

*"Examples of taxonomic gut microbiota composition. In the box are cited examples of bacteria belonging to Phyla Firmicutes and Bacteroidetes, representing 90% of gut microbiota" [4]*


## Литература
[1] Rebeca Martín, David Rios-Covian, Eugénie Huillet, Sandrine Auger, Sarah Khazaal, Luis G Bermúdez-Humarán, Harry Sokol, Jean-Marc Chatel, Philippe Langella, Faecalibacterium: a bacterial genus with promising human health applications, FEMS Microbiology Reviews, Volume 47, Issue 4, July 2023, fuad039, https://doi.org/10.1093/femsre/fuad039

[2] Eija Könönen, 250 - Anaerobic Cocci and Anaerobic Gram-Positive Nonsporulating Bacilli. 
Mandell, Douglas, and Bennett's Principles and Practice of Infectious Diseases (Eighth Edition),
Pages 2781-2786.e2. https://doi.org/10.1016/B978-1-4557-4801-3.00250-2

[3] R. Hutkins, Y.J. Goh, STREPTOCOCCUS | Streptococcus thermophilus, Encyclopedia of Food Microbiology (Second Edition), Academic Press, 2014, Pages 554-559, ISBN 9780123847331, https://doi.org/10.1016/B978-0-12-384730-0.00325-6

[4] Rinninella E, Raoul P, Cintoni M, Franceschi F, Miggiano GAD, Gasbarrini A, Mele MC. What is the Healthy Gut Microbiota Composition? A Changing Ecosystem across Age, Environment, Diet, and Diseases. Microorganisms. 2019 Jan 10;7(1):14. doi: 10.3390/microorganisms7010014. PMID: 30634578; PMCID: PMC6351938.



