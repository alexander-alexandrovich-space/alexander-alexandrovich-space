#include <ArduinoBLE.h>
#include <SCServo.h>

#define S_RXD 18
#define S_TXD 19

SCSCL servo;

void setup() {
  Serial.begin(115200);
  Serial1.begin(1000000, SERIAL_8N1, S_RXD, S_TXD);
  servo.pSerial = &Serial1;
  servo.WritePos(1,200,0,2000);
  servo.WritePos(2,340,0,2000);
  delay(1000);
  
  while (!Serial);

  BLE.begin();
  // start scanning for peripherals
  BLE.scanForUuid("19b10000-e8f2-537e-4f6c-d104768a1214");
}

void loop() {
  // check if a peripheral has been discovered
  BLEDevice peripheral = BLE.available();

  if (peripheral) {
    // discovered a peripheral, print out address, local name, and advertised service
    Serial.print("Found ");
    Serial.print(peripheral.address());
    Serial.print(" '");
    Serial.print(peripheral.localName());
    Serial.print("' ");
    Serial.print(peripheral.advertisedServiceUuid());
    Serial.println();


    // stop scanning
    BLE.stopScan();

    control(peripheral);

    // peripheral disconnected, start scanning again
    BLE.scanForUuid("19b10000-e8f2-537e-4f6c-d104768a1214");
  }
}

int servo_pos1=-1;
int servo_pos2=-2;
int accuracy_gor=5;
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
    //delay(1000);
    int read_pos1=-1, read_pos2=-1;
    int servo_speed1=1;

    for(servo_pos2=347 ; servo_pos2<=633; servo_pos2+=accuracy_vert){
      servo.WritePos(2,servo_pos2,0,2000);
      read_pos2=servo.ReadPos(2);
      
      while(servo_pos1<770){
        servo.WritePos(1, 777, 0, servo_speed1);
        servo_pos1=servo.ReadPos(1);
        servoCharacteristic1.writeValue((uint32_t)servo_pos1, 32);
        servoCharacteristic2.writeValue((uint32_t)read_pos2, 32);
        //Serial.print(servo_pos1);
        //Serial.print(",");
        //Serial.println(read_pos2);
        if(servo_pos1==-1){
        servo_pos1=203;
        servo_pos2=347;
        servo.WritePos(1,200,0,2000);
        servo.WritePos(2,340,0,2000);
        delay(1000);
        read_pos1=servo.ReadPos(1);
        read_pos2=servo.ReadPos(2);
        }
      }
      servo_pos2+=accuracy_vert;
      servo.WritePos(2,servo_pos2,0,2000);
      read_pos2=servo.ReadPos(2);
      while(servo_pos1>208){
        servo.WritePos(1, 203, 0, servo_speed1);
        servo_pos1=servo.ReadPos(1);
        servoCharacteristic1.writeValue((uint32_t)servo_pos1, 32);
        servoCharacteristic2.writeValue((uint32_t)read_pos2, 32);
        //Serial.print(servo_pos1);
        //Serial.print(",");
        //Serial.println(read_pos2);
        if(servo_pos1==-1){
        servo_pos1=203;
        servo_pos2=347;
        servo.WritePos(1,200,0,2000);
        servo.WritePos(2,340,0,2000);
        delay(1000);
        read_pos1=servo.ReadPos(1);
        read_pos2=servo.ReadPos(2);
        }
      }
      if(!peripheral.connected()){
        break;
      }
      
    }
    
  }

  Serial.println("Peripheral disconnected");
}