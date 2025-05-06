import processing.serial.*;
Serial myPort;  // The serial port
String myString = null;
float x, y, z;
int servo_pos1; // Азимут
int servo_pos2; // Угол места (элевейшн)
float distance;
ArrayList<PVector> pointCloud = new ArrayList<PVector>(); // Список для хранения точек

void setup() {
  // Список доступных последовательных портов
  printArray(Serial.list());
  // Открываем порт с нужной скоростью
  myPort = new Serial(this, Serial.list()[3], 115200);

  size(820, 820, P3D);
  noSmooth();
  background(0);
  stroke(255);
  strokeWeight(3);
}

void draw() {
  background(0); // Очищаем экран каждый кадр
  translate(width / 2, height / 2, 0); // Перемещаем начало координат в центр экрана
  rotateY(frameCount * 0.01); // Поворачиваем сцену для визуализации 3D эффекта

  // Отрисовываем все точки из облака
  stroke(255);
  strokeWeight(2);
  for (PVector p : pointCloud) {
    point(p.x, p.y, p.z);
  }

  // Читаем данные из последовательного порта
  while (myPort.available() > 0) {
    myString = myPort.readStringUntil(10); // Читаем до символа новой строки
    if (myString != null) {
      String[] q = splitTokens(myString, ",");
      if (q.length >= 3) {
        servo_pos1 = int(q[0]);  // Позиция сервопривода 1 (азимут)
        servo_pos2 = int(q[1]);  // Позиция сервопривода 2 (угол места)
        distance = float(q[2]);  // Дальность

        // Преобразование позиций сервоприводов в углы
        float azimuth_deg = ((servo_pos1 - 490) / 287.0) * 60.0;
        float elevation_deg = ((servo_pos2 - 490) / 143.0) * 15.0;

        // Преобразование углов в радианы
        float azimuth_rad = radians(azimuth_deg);
        float elevation_rad = radians(elevation_deg);

        // Преобразование сферических координат в декартовы
        x = distance * cos(elevation_rad) * sin(azimuth_rad);
        y = -distance * sin(elevation_rad);
        z = -distance * cos(elevation_rad) * cos(azimuth_rad);

        // Масштабирование координат для визуализации
        float scaleFactor = 0.5; // Измените по необходимости
        x *= scaleFactor;
        y *= scaleFactor;
        z *= scaleFactor;

        // Добавляем точку в облако точек
        pointCloud.add(new PVector(x, y, z));
      }
    }
  }
}
