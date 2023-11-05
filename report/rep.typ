

= Verschränkung

#let hbar = (sym.wj, move(dy: -0.08em, strike(offset: -0.55em, extent: -0.05em, sym.planck)), sym.wj).join()

Pauli Gleichung (ohne el. Potenzial):

$ i hbar diff_t phi =  1/(2m)[(p-q A)^2 - q hbar bold(sigma) dot B] phi $

Wobei:

- $p = -i hbar diff_z$
- $bold(sigma)$ ... Pauli Matrizen $sigma_1 = mat(0,1;1,0)$, $sigma_2 = mat(0,-i;i,0)$, $sigma_3 = mat(1,0;0,-1)$

Wir simulieren in 1D. Das Koordinatensystem wird so gwählt, dass B ausschliesslich in y Richtung zeigt. Damit lässt sich das magnetische Vektorpotenzial A recht einfach wählen als $ A_x (z) = integral_(z_min)^z B_y (z) dif x $

Es wird nur eine Komponente von $bold(sigma)$ benötigt, da $bold(sigma) dot B = sigma_2 B_y$.

Diese Gleichung wird nun diskretisiert:


$ i hbar diff_t phi &=  1/(2m)[(-i hbar diff_x)^2+ 2 (-i hbar diff_x) q A - (q A)^2 - g q hbar bold(sigma) dot B] phi \
 i hbar diff_t phi &=  1/(2m)[- hbar^2 diff_x^2 phi - 2i hbar q A diff_x phi - (q A)^2 phi - g q hbar bold(sigma)_2 B phi] \ 
 i hbar diff_t phi &=  1/(2m)[- hbar^2 (phi_(x-1)-2 phi_x + phi_(x+1))/(Delta x)^2 - 2 i hbar q A_x (phi_(x+1)-phi_(x-1))/(2 Delta x)- (q A_x)^2 phi_x - g q hbar bold(sigma)_2 B phi_x] \ 
 i hbar diff_t phi &= 
 1/(2m)[(-hbar^2/(Delta x)^2+i hbar q A_x / (Delta x))phi_(x-1) + 
 (-hbar^2/(Delta x)^2-i hbar q A_x / (Delta x) )phi_(x+1) + 
 (2 hbar^2/(Delta x)^2 - g q hbar sigma_2 B)] \ 
$
Dabei ist zu berücksichtigen dass $phi$ ein spinor ist also $phi = vec(phi arrow.t, phi arrow.b )$
Daher wurde die Gleichung wie folgt in Matrixschreibweise diskretisiert:
#pagebreak()
#set page("a4",flipped: true)



$ diff_t vec(dots.v, phi_(i-1) arrow.t,phi_(i-1) arrow.b,phi_(i) arrow.t,phi_(i) arrow.b,phi_(i+1) arrow.t,phi_(i+1) arrow.b,dots.v) = 1/(i hbar 2m)
mat(dots.down, dots.v, dots.v;
 dots,
 (-hbar^2/(Delta x)^2+i hbar q A_x/(Delta x)),0,
 (2 hbar^2/(Delta x)^2 - g q hbar sigma_2^[0,0] B),- g q hbar sigma_2^[0,1] B,
 (-hbar^2/(Delta x)^2-i hbar q A_x/(Delta x)), 0 ,  dots;
 dots,0,(-hbar^2/(Delta x)^2+i hbar q A_x/(Delta x)),
  - g q hbar sigma_2^[1,0] B, 2 hbar^2/(Delta x)^2 - g q hbar sigma_2^[1,1] B,
   0, (-hbar^2/(Delta x)^2-i hbar q A_x /(Delta x))
      )  
  vec(dots.v, phi_(i-1) arrow.t,phi_(i-1) arrow.b,phi_(i) arrow.t,phi_(i) arrow.b,phi_(i+1) arrow.t,phi_(i+1) arrow.b, dots.v) $ 
#pagebreak()
#set page("a4",flipped: false)
Das Vektorpotenzial wurde wie folgt gewählt (B0 = 1T) mit linearer Übergangsphase:
#figure(image("potential.svg", width: 70%))

Für das Wellenpaket wurde ein Gaussches Wellenpaket mit spin up gewählt:
$ phi =  vec((1/(2 pi a^2))^(1/4) exp(-(z-z_0)^2/(4 a^2))  exp(i k_0 z),0) $ 

= Wellenfunktion für Zwei Partikel

* Ab hier bin ich sehr unsicher ob das so stimmt, ich habe mir ein paar dinge zu Zwei-Partikel Schrödinger-Gleichungen angesehen und dies versucht auf die Pauli-Gleichung zu übertragen *


Als nächstes habe ich nun probiert die Wellenfunktion für zwei Teilchen zu simulieren. Dazu hätte ich einen - Zweiteilchen Spinor angesetzt: $phi(x_1,x_2)$ und für die Pauli gleichung einfach die Energieterme der beiden Partikel addiert:
$ i hbar diff_t  phi(x_1,x_2)  =  1/(2m)[(-i hbar diff_x_1)^2 + (-i hbar diff_x_2)^2+ \ 
(-i hbar diff_x_1) q A[x_1] +  (-i hbar diff_x_2) q A[x_2] \ 
- (q A[x_1])^2 - (q A[x_2])^2 \
- g q hbar bold(sigma) dot B[x_1] - g q hbar bold(sigma) dot B[x_2] ] phi(x_1,x_2) 
$

Diese Gleichung habe ich in ähnlicher Weise diskretisiert (eben 2-dimensional).
Allerdings bin ich mir hier nicht ganz sicher wie ich weiter fortfahren sollte.
Ich kann "einzelne" Teilchen einsetzen (die eine konstante Funktion in einer Komponente sind (mit passender Skalierung sodass die Norm 1 wird))
und erhalte dabei das gleiche Ergebnis wie bei zwei Teilchen.
Ich bin mir aber nicht ganz sicher, wie das mit der verschränkung funktioniert.
Wenn ich einen "reinen" up zustand mit einem "reinen" down zustand multipliziere kommt insgesamt 0 heraus (weil die andere spin Komponente eben 0 ist).

Um die Lösung für ein Partikel zu erhalten, habe ich über die andere (absolut-quadrierte) Komponente integriert:
$ |phi(x_1)|^2 = integral |phi(x_1,x_2)|^2 dif x_2$