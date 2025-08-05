# 🚀 alexander-alexandrovich-space

Добро пожаловать в **alexander-alexandrovich-space** — сборник скриптов, моделей и отчётов по различным  космическихим системам, электронике и цифровой обработки сигналов.

---

## 🌐 Simulation of Intersatellite Network

В этой части репозитория собраны MATLAB-скрипты для моделирования и анализа межспутниковой сети:

- **simulation of intersatellite network.m** и **glonass.m**  
  - Отрисовывают орбиты и размещают спутники на орбитах заданной конфигурации (их количество, высотный диапазон, количество и плоскости можно регулировать)  
  - Вычисляют истинные аномалии каждого спутника.  
  - Решают задачи организации связи для заданного момента времени по заданным параметрам приемо-передающих устройств.
  - Есть модель протокола OLSR для меш-сети.
  - Есть алгоритм поиска изолированных узлов и кластеров спутников.
  - Есть модель непосредственной модуляции сигнала при помощи CSS, используемой в технологии LoRa.
  - Рассчитывается теоретическая дальность связи Создается матрица случайного канального шума от 0 до 5 дБ для оценки дальности.
  - Строится матрица коммутаций, задержек распространения и граф сети на основе пересчитанной реальной дальности с уже внесенным канальным шумом. 
  - Поиск наиболее длинного маршрута в сети, который называется наихудшим, поскольку в этой модели пока нет симуляции загрузки сети и каналы оцениваются только по отношению сигнал-шум, это будет показательная оценка, так как более длинного маршрута в незагруженной сети быть не может.
  - Есть моделирование передачи с учетом отношения сигнал-шум и расчет задержек передачи на основе типовых значений скорости передачи и обработки данных меш-узлами.

  - Для корректной работы необходим файл:
    - `earthshmap.jpg` — текстура поверхности Земли.
  - Постер с конференции **SPBOPEN 2025** в качестве инфографики.
  - Работа отмечена грамотой за лучший доклад на международной конференции "Ломоносов" в 2025 году
---

## 🛠 WARD

Папка `WARD` содержит наработки по различным механизмам и 3D-моделям:

- Скрипты для управления и визуализации:
  - Моделирование работы 3D-лазерного лидара из двух микроконтроллеров.
  - Пост-процессинг и отрисовка облаков точек в Processing.
- **3D-модель корпуса** в формате `.stl` для дальнейшей интеграции в CAD.

---

## 🔌 Power Supply

В папке `power_supply` — **курсовая работа** по разработке импульсного источника питания:

1. Исходные схемы и модели SPICE.
2. Отчёт с анализом характеристик:
   - КПД при разных нагрузках.
   - Анализ пульсаций.
---

## 📊 DSP

Папка `DSP` содержит отчёт по лабораторным работам курса **Цифровой обработки сигналов**:

- Анализ фильтров (IIR/FIR).
- Спектральный анализ сигналов.
- Преобразование Гильберта

---

## 💾 Курсовая работа

Архив `Курсовая работа` включает курсовую работу по **цифровым устройствам и микропроцессорам**, в которой реализована игра «Больше — Меньше» для микроконтроллера **AT32**:

- **Генерация случайного числа**  
  Считывание напряжения с термодатчика и преобразование в целое число от 0 до 100.
- **Ввод с 16-символьной клавиатуры**  
  Пользователь по UART вводит угадываемое число.
- **Логика игры**  
  После каждой попытки микроконтроллер сообщает, больше или меньше введённое число по сравнению с загаданным.
- **Конец игры**  
  При угадывании выводится поздравление и количество затраченных попыток.
---
# 💻 Tech Stack:
![C](https://img.shields.io/badge/c-%2300599C.svg?style=flat&logo=c&logoColor=white) ![C++](https://img.shields.io/badge/c++-%2300599C.svg?style=flat&logo=c%2B%2B&logoColor=white) ![Bash Script](https://img.shields.io/badge/bash_script-%23121011.svg?style=flat&logo=gnu-bash&logoColor=white) ![Windows Terminal](https://img.shields.io/badge/Windows%20Terminal-%234D4D4D.svg?style=flat&logo=windows-terminal&logoColor=white) ![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=flat&logo=github&logoColor=white) ![Git](https://img.shields.io/badge/git-%23F05033.svg?style=flat&logo=git&logoColor=white) ![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=flat&logo=docker&logoColor=white) ![Raspberry Pi](https://img.shields.io/badge/-Raspberry_Pi-C51A4A?style=flat&logo=Raspberry-Pi) ![Ansible](https://img.shields.io/badge/ansible-%231A1918.svg?style=flat&logo=ansible&logoColor=white) ![Arduino](https://img.shields.io/badge/-Arduino-00979D?style=flat&logo=Arduino&logoColor=white)
# 📊 GitHub Stats:
![](https://github-readme-stats.vercel.app/api/top-langs/?username=alexander-alexandrovich-space&theme=dark&hide_border=false&include_all_commits=false&count_private=false&layout=compact)

### 🔝 Top Contributed Repo
![](https://github-contributor-stats.vercel.app/api?username=alexander-alexandrovich-space&limit=5&theme=dark&combine_all_yearly_contributions=true)

---

> **Автор:** Александр Александрович Попов 
> **Контакты:** alexander_popov_work@mail.ru
