#include <ArduinoBLE.h>
#include <SCServo.h>
#include <math.h>

#define S_RXD 18
#define S_TXD 19

SCSCL servo;

 int servo_pos3=-3;
 int servo_pos4=-4;

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
  servo.WritePos(1, 512, 0, 0);
  servo.WritePos(2, 512, 0, 0);
  servo.WritePos(3, 512, 0, 0);
  servo.WritePos(4,512,0,0);
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


  BLECharacteristic servoCharacteristic3 = peripheral.characteristic("19b10001-e8f2-537e-4f6c-d104768a1204");
  BLECharacteristic servoCharacteristic4 = peripheral.characteristic("19b10001-e8f2-537e-4f6c-d104768a1205");

    unsigned long laserVal;

  while (peripheral.connected()) {
    servo_pos3=servo.ReadPos(3);
    servo_pos4=servo.ReadPos(4);
    uint32_t rnd = random(0,1000);
    servoCharacteristic3.writeValue(rnd);
    servoCharacteristic4.writeValue(rnd);

    laserCharacteristic.readValue((uint8_t*)&laserVal, sizeof(unsigned long));   // читаем 4 байта
    Serial.println(laserVal);
    servo.WritePos(3, 512+(atan(167/(laserVal))*(180/PI))*(800/180) , 0, 0); // 3 - ЛЕВЫЙ, ЕСЛИ СМОТРЕТЬ В СТОРОНУ ЛАЗЕРА
    servo.WritePos(4, 512-(atan(167/(laserVal))*(180/PI))*(800/180) , 0, 0); // 4 - ЛЕВЫЙ, ЕСЛИ СМОТРЕТЬ В СТОРОНУ ЛАЗЕРА
    }
    

  Serial.println("Peripheral disconnected");
}