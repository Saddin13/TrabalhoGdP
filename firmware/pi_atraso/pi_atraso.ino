// PI discreto com atraso simulado (buffer circular)
// Planta: RC de 2a ordem (R1=10k, C1=10uF, R2=100k, C2=1uF)
// Saida PWM: D9   |   Leitura: A0
// Usado no Tinkercad Circuits (Arduino Uno R3)

const float T   = 0.01;   // periodo de amostragem: 10 ms
const float Kp  = 2.0;    // ganho proporcional
const float Ki  = 1.0;    // ganho integral
const float REF = 2.5;    // referencia (V)
const int   N   = 0;      // ATRASO em amostras: tau = N*T (critico observado: N = 15)

const int NMAX = 250;     // atraso maximo: 250 amostras = 2,5 s
uint8_t buf[NMAX + 1];    // buffer circular com as acoes de controle
int idx = 0;
float I = 0;              // termo integral
unsigned long t0;
int cont = 0;

void setup() {
  pinMode(9, OUTPUT);
  Serial.begin(9600);
  for (int i = 0; i <= NMAX; i++) buf[i] = 0;
  t0 = millis();
}

void loop() {
  // espera completar o periodo de amostragem
  if (millis() - t0 < (unsigned long)(T * 1000)) return;
  t0 += (unsigned long)(T * 1000);

  // leitura da saida da planta (0 a 5 V)
  float y = analogRead(A0) * 5.0 / 1023.0;

  // controlador PI
  float e = REF - y;
  I += Ki * e * T;
  I = constrain(I, -5, 5);                 // anti-windup
  float u = constrain(Kp * e + I, 0, 5);

  // atraso: guarda a acao atual e aplica a de N amostras atras
  buf[idx] = (uint8_t)(u / 5.0 * 255);
  int saida = (idx - N + NMAX + 1) % (NMAX + 1);
  analogWrite(9, buf[saida]);
  idx = (idx + 1) % (NMAX + 1);

  // envia a cada 100 ms: referencia, saida, controle
  if (++cont >= 10) {
    cont = 0;
    Serial.print(REF); Serial.print(" ");
    Serial.print(y);   Serial.print(" ");
    Serial.println(u);
  }
}
