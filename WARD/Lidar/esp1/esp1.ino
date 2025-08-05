#include <ArduinoBLE.h>
#include <SCServo.h>
#include <math.h>

#define S_RXD 18
#define S_TXD 19

SCSCL servo;

void setup() {
  Serial.begin(115200);
  Serial1.begin(1000000, SERIAL_8N1, S_RXD, S_TXD);
  servo.pSerial = &Serial1;
  
  while (!Serial);

  BLE.begin();
  // start scanning for peripherals
  BLE.scanForUuid("19b10000-e8f2-537e-4f6c-d104768a1214");
}

void loop() {
  servo.WritePos(3, 512, 0, 0);
  servo.WritePos(4,512,0,0);
  delay(1000);
  // check if a peripheral has been discovered
  BLEDevice peripheral = BLE.available();

  if (peripheral) {
    // discovered a peripheral, print out address, local name, and advertised service
    Serial.print("Found ");
    Serial.print(peripheral.address());
    Serial.print(" ");
    Serial.print(peripheral.localName());
    Serial.print(" ");
    Serial.print(peripheral.advertisedServiceUuid());
    Serial.println();


    // stop scanning
    BLE.stopScan();

    control(peripheral);

    // peripheral disconnected, start scanning again
    //BLE.scanForUuid("19b10000-e8f2-537e-4f6c-d104768a1214");
    BLE.scanForUuid("19b10000-e8f2-537e-4f6c-d104768a1214");

  }
}

int servo_pos1=-1;
int servo_pos2=-2;
//int accuracy_gor=10;
int accuracy_vert=5;


void control(BLEDevice peripheral) {
  // connect to the peripheral
  Serial.println("Connecting ...");

  if (peripheral.connect()) {
    Serial.println("Connected");
  } else {
    Serial.println("Failed to connect!");
    return;
  }

  // discover peripheral attributes
  Serial.println("Discovering attributes ...");
  if (peripheral.discoverAttributes()) {
    Serial.println("Attributes discovered");
  } else {
    Serial.println("Attribute discovery failed!");
    peripheral.disconnect();
    return;
  }

  // retrieve the LED characteristic
  BLECharacteristic servoCharacteristic1 = peripheral.characteristic("19b10001-e8f2-537e-4f6c-d104768a1201");
  BLECharacteristic servoCharacteristic2 = peripheral.characteristic("19b10001-e8f2-537e-4f6c-d104768a1202");
  BLECharacteristic laserCharacteristic = peripheral.characteristic("19B10001-E8F2-537E-4F6C-D104768A1203");

  if (!servoCharacteristic1 && !servoCharacteristic2) {
    Serial.println("Peripheral does not have servo characteristic!");
    peripheral.disconnect();
    return;
  } else if (!servoCharacteristic1.canWrite() && !servoCharacteristic2.canWrite()) {
    Serial.println("Peripheral does not have a writable servo characteristic!");
    peripheral.disconnect();
    return;
  }

  while (peripheral.connected()) {
    int read_pos1=-1, read_pos2=-1;
    int servo_speed1=1;

    servo.WritePos(1,512,0,0);
    servo.WritePos(2,512,0,0);
    delay(1000);
    servo_pos1=servo.ReadPos(1);
    servo_pos2=servo.ReadPos(2);
    servoCharacteristic1.writeValue((uint32_t)servo_pos1, 32);
    servoCharacteristic2.writeValue((uint32_t)servo_pos2, 32);

    unsigned long laserVal;

    laserCharacteristic.readValue((uint8_t*)&laserVal, sizeof(unsigned long));   // читаем 4 байта
    Serial.println(laserVal);
    servo.WritePos(3, 512+(atan(167/laserVal)*(180/PI))*(800/180) , 0, 0); // 3 - ЛЕВЫЙ, ЕСЛИ СМОТРЕТЬ В СТОРОНУ ЛАЗЕРА
    servo.WritePos(4, 512-(atan(167/laserVal)*(180/PI))*(800/180) , 0, 0); // 4 - ЛЕВЫЙ, ЕСЛИ СМОТРЕТЬ В СТОРОНУ ЛАЗЕРА

     servo.WritePos(1,448,0,0);
     servo.WritePos(2,467,0,0);

    delay(1000);

    servo_pos1=servo.ReadPos(1);
    servo_pos2=servo.ReadPos(2);

    int servo_pos_nach_1 = servo_pos1;
    int servo_pos_nach_2 = servo_pos2;

    for(servo_pos2=467 ; servo_pos2<=557; servo_pos2+=accuracy_vert){
      servo.WritePos(2,servo_pos2,0,0);
      delay(100);
      read_pos2=servo.ReadPos(2);
      
      while(servo_pos1<576){
        servo.WritePos(1, 586, 0, servo_speed1);
        servo_pos1=servo.ReadPos(1);
        servoCharacteristic1.writeValue((uint32_t)servo_pos1, 32);
        servoCharacteristic2.writeValue((uint32_t)read_pos2, 32);
        laserCharacteristic.readValue((uint8_t*)&laserVal, sizeof(unsigned long));   // читаем 4 байта

        Serial.print("Координата x в дел: ");
        Serial.print(servo_pos1);
        Serial.print(",");
        Serial.print("Координата y в дел: ");
        Serial.print(servo_pos2);
        Serial.print(",");
        Serial.print("Дальность в мм: ");
        Serial.print(laserVal);
        Serial.println(",");

        //Serial.println(servo.ReadPos(3));
        //Serial.println(servo.ReadPos(4));
        if(servo_pos1==-1){
        servo_pos1=448;
        servo_pos2=467;
        servo.WritePos(1,448,0,0);
        servo.WritePos(2,467,0,0);
        delay(1000);
        read_pos1=servo.ReadPos(1);
        read_pos2=servo.ReadPos(2);
        }
        laserCharacteristic.readValue((uint8_t*)&laserVal, sizeof(unsigned long));   // читаем 4 байта
        //Serial.println(laserVal);
      }
      servo_pos2+=accuracy_vert;
      servo.WritePos(2,servo_pos2,0,0);
      read_pos2=servo.ReadPos(2);
      while(servo_pos1>448){
        servo.WritePos(1, 438, 0, servo_speed1);
        servo_pos1=servo.ReadPos(1);
        servoCharacteristic1.writeValue((uint32_t)servo_pos1, 32);
        servoCharacteristic2.writeValue((uint32_t)read_pos2, 32);

        Serial.print("Координата x в дел: ");
        Serial.print(servo_pos1);
        Serial.print(",");
        Serial.print("Координата y в дел: ");
        Serial.print(servo_pos2);
        Serial.print(",");
        Serial.print("Дальность в мм: ");
        Serial.print(laserVal);
        Serial.println(",");

       // Serial.println(servo.ReadPos(3));
        //Serial.println(servo.ReadPos(4));
        if(servo_pos1==-1){
        servo_pos1=448;
        servo_pos2=467;
        servo.WritePos(1,448,0,0);
        servo.WritePos(2,467,0,0);
        delay(1000);
        read_pos1=servo.ReadPos(1);
        read_pos2=servo.ReadPos(2);
        

        }
        //laserCharacteristic.readValue((uint8_t*)&laserVal, sizeof(unsigned long));   // читаем 4 байта
        //Serial.println(laserVal);
      }
      if(!peripheral.connected()){
        break;
      }
      
    } // конец фора
    Serial.println("");
    Serial.print("макс Размер 1 в мм: ");
    Serial.print(tan(((servo_pos1-servo_pos_nach_1)*(PI/800))/2)*1200);
    Serial.print(",");
    Serial.print("макс Размер 2 в мм: ");
    Serial.print(tan(((servo_pos2-servo_pos_nach_2)*(PI/800))/2)*1200);
    Serial.println(",");

    delay(2000);
    
  }

  Serial.println("Peripheral disconnected");
}