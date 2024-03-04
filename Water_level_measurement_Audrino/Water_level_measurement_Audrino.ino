// #define ledPin 6
#define sensorPin A5
unsigned long CurrentTime;
unsigned long StartTime;
float time_in_sec;
bool start_clock;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(1200);
  start_clock = true;
  Serial.println("Time [s], liquid_level_reading");
  
}

void loop() {
  // put your main code here, to run repeatedly:

  if(start_clock) 
  {
    StartTime = millis();
    CurrentTime = StartTime;
    start_clock = false;
  }
  else CurrentTime = millis();
  time_in_sec = (CurrentTime - StartTime)/1000.0;
  Serial.print(time_in_sec);
  Serial.print(",");
  int sensorValue = analogRead(sensorPin);
  Serial.println(sensorValue);
  // Serial.print("\n");
 
}
// 5v
// empty - 508
// 1 inch - 514
// 2 inch - 533
// 3 inch - 554
// 4 inch - 576
// 5 inch - 606
// 6 inch - 623
// 7 inch - 649
// 8 inch - 682
// 9 inch - 704
// 10 inch - 754
// 11 inch - 805 

