int prescaler = 256; // set this to match whatever prescaler value you set in CS registers below
int speed = 200; // average the pressure measurements from this many samples

// intialize values for the PWM duty cycle set by pots
float potDC1 = 0;
float potDC2 = 0;
float potDC3 = 0;
float potDC4 = 0;

byte scaling = 0xFF;

bool automode = false;

void setup() {

  Serial.begin(9600);

  // input pins for valve switches
  pinMode(50, INPUT);
  pinMode(51, INPUT);
  pinMode(52, INPUT);
  pinMode(53, INPUT);

  // output pins for valve PWM
  pinMode(5, OUTPUT);
  pinMode(6, OUTPUT);
  pinMode(7, OUTPUT);
  pinMode(8, OUTPUT);

  int eightOnes = 255;  // this is 11111111 in binary
  TCCR3A &= ~eightOnes;   // this operation (AND plus NOT), set the eight bits in TCCR registers to 0 
  TCCR3B &= ~eightOnes;
  TCCR4A &= ~eightOnes;
  TCCR4B &= ~eightOnes;

  // set waveform generation to frequency and phase correct, non-inverting PWM output
  TCCR3A = _BV(COM3A1);
  TCCR3B = _BV(WGM33) | _BV(CS32);
  
  TCCR4A = _BV(COM4A1) | _BV(COM4B1) | _BV(COM4C1);
  TCCR4B = _BV(WGM43) | _BV(CS42);
}

void pPWM(float pwmfreq, float pwmDC1, float pwmDC2, float pwmDC3, float pwmDC4) {

  // set PWM frequency by adjusting ICR (top of triangle waveform)
  ICR3 = F_CPU / (prescaler * pwmfreq * 2);
  ICR4 = F_CPU / (prescaler * pwmfreq * 2);
  
  // set duty cycles
  OCR3A = (ICR4) * (pwmDC1 * 0.01);
  OCR4A = (ICR4) * (pwmDC2 * 0.01);
  OCR4B = (ICR4) * (pwmDC3 * 0.01);
  OCR4C = (ICR4) * (pwmDC4 * 0.01);
}

void loop() {

  potDC1 = 0; potDC2 = 0; potDC3 = 0; potDC4 = 0;

  if (Serial.available() > 0) {                    
    scaling = Serial.read();   
    automode = true; 
  }

  // if statement for manual switch override
  if (digitalRead(50) == HIGH) { potDC1 = analogRead(A1)*100.0/1024.0 * ((scaling & 0x03) / 3.0); }
  if (digitalRead(51) == HIGH) { potDC2 = analogRead(A2)*100.0/1024.0 * (((scaling & 0x0C) >> 2) / 3.0); }

  float tmp = analogRead(A3)*100.0/1024.0;
  if (digitalRead(52) == HIGH) { potDC3 = tmp * (((scaling & 0x30) >> 4) / 3.0); }  
  if (digitalRead(53) == HIGH) { potDC4 = tmp * (((scaling & 0xC0) >> 6) / 3.0); }  
  
  float PWMfq = analogRead(A4)*100.0/1024.0; // scale values from pot to 0 to 100, which gets used for frequency (Hz)
  PWMfq = round(PWMfq/5)*5+1; //1 to 91 Hz in increments of 5 (rounding helps to deal with noisy pot)

  if (automode) {
    PWMfq = 56;
  }

  // update PWM output based on the above values from pots
  pPWM(PWMfq,potDC1,potDC2,potDC3,potDC4);

  // transfer function for sensor Honeywell ASDXRRX100PGAA5 (100 psi, 5V, A-calibration)
  // Vout = 0.8*Vsupply/(Pmax - Pmin)*(Papplied - Pmin) + 0.1*Vsupply
  // Rearrange to get: Papplied = (Vout/Vsupply - 0.1)*(Pmax - Pmin)/0.8 + Pmin;

  // to convert to pressure reading into PSI
  // float p_psi = (value/1024.0 - 0.1)*100.0/0.8;  
  // to convert to pressure reading into kPa
  // float p_kPa = (value/1024.0 - 0.1)*100.0/0.8 * 6.89475729;  

  // print pressure readings
  long c1 = 0,c2 = 0,c3 = 0,c4 = 0;
  
  for(int i = 0;i < speed;i++)
  {
    c1 = c1 + analogRead(A8);
    c2 = c2 + analogRead(A9);
    c3 = c3 + analogRead(A10);
    c4 = c4 + analogRead(A11);    
  }
  Serial.print(c1 / speed); Serial.print("\t");
  Serial.print(c2 / speed); Serial.print("\t");
  Serial.print(c3 / speed); Serial.print("\t");
  Serial.print(c4 / speed); Serial.print("\t");
  Serial.print(PWMfq); Serial.print("\n");
}
