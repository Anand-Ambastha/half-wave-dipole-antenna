# Theory of the Half-Wave Dipole

This note covers the theory used in the project: the current on the wire, the radiated fields, impedance, S11 and VSWR, and how the CST and Python simulations work.

Symbols: f is frequency, λ = c/f is wavelength, k = 2π/λ, η ≈ 377 Ω is the impedance of free space, L is the total length of the dipole and a is the wire radius.

## 1. The half-wave dipole

A dipole is a straight wire fed at the centre across a small gap. When its total length is close to half a wavelength,

$$L \approx \frac{\lambda}{2},$$

the current on the wire forms a standing wave and the input reactance is close to zero. The antenna is then resonant and takes in power easily. Its polarisation is linear, with the electric field parallel to the wire.

At 2.4 GHz, λ = 124.9 mm, so λ/2 = 62.4 mm. The physical resonant length is a little shorter, as explained in section 5.

## 2. Current distribution

A dipole can be pictured as an open-ended transmission line whose two wires have been opened out. The current is then a standing wave that drops to zero at both tips:

$$I(z) = I_0 \sin\!\left[k\left(\frac{L}{2} - |z|\right)\right], \qquad |z| \le \frac{L}{2}.$$

For L = λ/2 this reduces to I(z) = I₀ cos(kz): largest at the feed point and zero at the ends. This is a good approximation for thin wires (a much smaller than λ) and is the basis of the closed-form results below. The Python program does not assume this shape. It calculates the current numerically, and `results/current.png` shows the result.

## 3. Radiated field and pattern

Taking the far-field part of the field (the part that falls off as 1/r), the electric field of a centre-fed dipole is

$$E_\theta = j\,\frac{\eta I_0 e^{-jkr}}{2\pi r}\;
\frac{\cos\!\left(\frac{kL}{2}\cos\theta\right) - \cos\!\left(\frac{kL}{2}\right)}{\sin\theta},
\qquad H_\phi = \frac{E_\theta}{\eta},$$

where θ is measured from the wire axis. For L = λ/2 this becomes

$$E_\theta = j\,\frac{\eta I_0 e^{-jkr}}{2\pi r}\;\frac{\cos\!\left(\frac{\pi}{2}\cos\theta\right)}{\sin\theta}.$$

The radiation intensity (power per unit solid angle) is

$$U(\theta) = \frac{\eta I_0^2}{8\pi^2}\left[\frac{\cos\!\left(\frac{\pi}{2}\cos\theta\right)}{\sin\theta}\right]^2.$$

What this means for the pattern:

- In the E-plane (any plane that contains the wire) the pattern is a figure of eight. It is largest at θ = 90° and has nulls along the wire.
- In the H-plane (the plane perpendicular to the wire) the field does not depend on the angle φ, so the pattern is a circle.
- In 3D the pattern is a doughnut around the wire.
- The half-power beamwidth in the E-plane is about 78°. A very short dipole has 90°.

## 4. Radiated power, radiation resistance and directivity

The total radiated power is

$$P_{rad} = \int_0^{2\pi}\!\!\int_0^{\pi} U \sin\theta\, d\theta\, d\phi
= \frac{\eta I_0^2}{4\pi}\int_0^{\pi}\frac{\cos^2\!\left(\frac{\pi}{2}\cos\theta\right)}{\sin\theta}\,d\theta
\approx 36.5\, I_0^2 \ \text{W}.$$

The integral equals half of Cin(2π), which is about 1.22, where Cin is the modified cosine integral.

The radiation resistance, referred to the current maximum (which is also the feed current here), is

$$R_r = \frac{2P_{rad}}{I_0^2} \approx 73.1\ \Omega.$$

The directivity is

$$D_0 = \frac{4\pi U_{max}}{P_{rad}} \approx 1.643 \quad (2.15\ \text{dBi}).$$

That is about 0.4 dB more than a very short dipole (1.5, or 1.76 dBi) and 2.15 dB more than an isotropic radiator.

Losses in the wire reduce the gain:

$$G = e_{cd}\,D_0, \qquad e_{cd} = \frac{R_r}{R_r + R_{loss}}.$$

The skin-effect loss resistance of a round wire, referred to the feed current, is

$$R_{loss} = \frac{R_s}{2\pi a}\cdot\frac{1}{|I_{feed}|^2}\int_{-L/2}^{L/2} |I(z)|^2\,dz,
\qquad R_s = \sqrt{\frac{\pi f \mu_0}{\sigma}}.$$

`dipole_sim.py` evaluates this from the calculated current. For copper (σ = 5.8 × 10⁷ S/m) at 2.4 GHz the loss is very small and the efficiency is about 99.9 %. With a perfect conductor (PEC) it is 100 %.

## 5. Input impedance and the shortening of the resonant length

For a thin dipole of exactly λ/2 the induced-EMF method gives

$$Z_{in} \approx 73.1 + j\,42.5\ \Omega.$$

The positive reactance means the wire is a little too long to be resonant. Making it shorter lowers the reactance. For a thin wire the reactance is zero at

$$L \approx 0.47 \text{ to } 0.48\,\lambda,$$

where the resistance is slightly below 73 Ω. A common rule of thumb is L ≈ 0.95 × λ/2. The shortening comes from the end effect (the current does not fall to zero at the tip exactly as the ideal line model assumes) and from the finite thickness of the wire.

Wire thickness is often described by

$$\Omega = 2\ln\!\left(\frac{L}{a}\right).$$

A thicker wire has a lower characteristic impedance as a transmission line. The reactance then changes more slowly with frequency, so the bandwidth is larger, and the resonant length becomes shorter. The radius sweep shows both effects: with the Python program the resonant length goes from 58.8 mm at a = 0.25 mm to 57.1 mm at a = 1 mm.

The feed gap also matters. A real gap adds a little capacitance and moves the resonance slightly, which is why the gap width is included in the parameter study.

## 6. S11, VSWR and bandwidth

For a reference impedance Z₀ = 50 Ω the reflection coefficient is

$$\Gamma = S_{11} = \frac{Z_{in} - Z_0}{Z_{in} + Z_0},$$

$$S_{11}\,[\text{dB}] = 20\log_{10}|\Gamma|, \qquad
\mathrm{VSWR} = \frac{1 + |\Gamma|}{1 - |\Gamma|}.$$

The fraction of power accepted by the antenna is 1 − |Γ|². A common design target is S11 below −10 dB, which means VSWR below 2 and at least 90 % of the power accepted.

Example: Z_in = 75.5 Ω gives |Γ| = 0.2, S11 = −14 dB and VSWR = 1.5. A 73 Ω dipole on a 50 Ω line is therefore already reasonably matched without a matching network.

The fractional bandwidth is the frequency range in which S11 stays below −10 dB:

$$\mathrm{BW}\,[\%] = \frac{f_{high} - f_{low}}{f_0}\times 100.$$

A dipole is a balanced antenna and coaxial cable is unbalanced, so a real build would use a balun to stop current flowing on the outside of the cable. The simulations here drive the dipole directly with an ideal port and do not model a balun.

## 7. Far-field distance

Patterns are defined in the far field, beyond the distance

$$r > \frac{2D^2}{\lambda}, \qquad D = L.$$

For this antenna that is about 62 mm (roughly 0.5 λ) at 2.4 GHz. Both simulations calculate the far field directly from the currents, so the pattern does not depend on the observation distance.

## 8. Simulation methods

### 8.1 CST Microwave Studio

The time-domain solver in CST uses the Finite Integration Technique, which is closely related to FDTD, on a hexahedral mesh.

1. The two PEC cylinders (2 mm gap between them) sit inside a computational box with open boundaries, about λ/4 or more away from the antenna.
2. A discrete port across the gap sends in a broadband pulse in a 50 Ω system.
3. Maxwell's equations are stepped in time until the energy in the box has decayed. The accuracy setting is −40 dB.
4. S-parameters and impedance come from a Fourier transform of the time signals, so one run covers the whole 1.5 to 3.5 GHz range.
5. Far-field monitors at 2.2, 2.4 and 2.6 GHz give the pattern, directivity and gain.

The main accuracy control is the mesh near the thin wire and the gap. The macro asks for 20 cells per wavelength, and the mesh should be checked and refined around the gap.

### 8.2 Python method of moments

`python/dipole_sim.py` solves the electric field integral equation for a thin straight wire. The tangential electric field on the wire surface must be zero except at the feed:

$$E_z^{s}(z) + E_z^{i}(z) = 0.$$

The scattered field comes from the unknown current I(z):

$$E_z^{s} = -j\omega A_z - \frac{\partial \Phi}{\partial z}, \qquad
A_z = \mu_0\!\int I(z')\,G\,dz', \qquad
\Phi = -\frac{1}{j\omega\varepsilon_0}\int \frac{dI}{dz'}\,G\,dz',$$

with the thin-wire kernel

$$G(z,z') = \frac{e^{-jkR}}{4\pi R}, \qquad R = \sqrt{(z-z')^2 + a^2}.$$

The steps are:

1. Divide the wire into N segments and write the current as a sum of triangular basis functions T_n(z). The current is zero at both ends by construction.
2. Test the equation with the same triangles (the Galerkin method). This gives a matrix equation

   $$\mathbf{Z}\,\mathbf{I} = \mathbf{V}, \qquad
   Z_{mn} = j\omega\mu_0\langle T_m, G\,T_n\rangle
   + \frac{1}{j\omega\varepsilon_0}\langle T_m', G\,T_n'\rangle .$$

3. The excitation V_m is the integral of T_m times the applied field. A 1 V source is spread uniformly over the 2 mm gap, which matches the CST discrete port.
4. Solve for the currents and take Z_in = 1 V / I_feed.
5. S11 and VSWR follow from Z_in. The far field comes from integrating the current against the phase factor e^{jkz′cosθ}. Directivity comes from integrating the radiation intensity over the sphere. The beamwidth is read from the pattern, and the ohmic loss from the skin resistance.

Checks that were done: for a very thin wire the program gives 73 + j42 Ω for a λ/2 dipole, the radiated power calculated from the current equals the power delivered at the feed, and the reactance of a short dipole is within about 8 % of theory.

Limits: the thin-wire kernel assumes the radius is much smaller than the wavelength and that segments are not shorter than the radius. For the 1 mm wire the impedance away from resonance can be off by 10 to 25 % and changes by a few percent with mesh density, so the CST results should be used as the reference. The pattern, directivity and trend studies are not sensitive to this.

## 9. Expected trends in the parameter studies

| Parameter changed | Expected effect |
|---|---|
| Length L | Longer wire gives a lower resonant frequency. The reactance changes sign at resonance and the input resistance rises with length. |
| Wire radius a | Larger radius gives more bandwidth and a shorter resonant length. |
| Feed gap g | Larger gap gives a small shift in resonance and impedance. |
| Material | Copper instead of PEC lowers the efficiency by only about 0.1 %. |
| Frequency | The reactance is capacitive below resonance and inductive above it. |

## 10. Comparison of theory and simulation

| Quantity | Theory (thin half-wave dipole) | Python MoM (a = 1 mm) | CST |
|---|---|---|---|
| Input resistance | 73.1 Ω | 75.5 Ω at resonance | to be filled |
| Resonant length | about 0.47 to 0.48 λ | 0.457 λ (57.1 mm) | to be filled |
| Directivity | 2.15 dBi | 2.14 dBi | to be filled |
| E-plane beamwidth | 78° | 78.5° | to be filled |
| Efficiency (PEC) | 100 % | 100 % | to be filled |
| S11 at 2.4 GHz | not applicable | −13.8 dB | to be filled |

The theory column assumes a very thin, perfect conductor. The 1 mm wire used here has a/λ = 0.008, which is fairly thick, and that explains the shorter resonant length.

## References

1. C. A. Balanis, *Antenna Theory: Analysis and Design*, Wiley.
2. W. L. Stutzman and G. A. Thiele, *Antenna Theory and Design*, Wiley.
3. J. D. Kraus and R. J. Marhefka, *Antennas for All Applications*, McGraw-Hill.
4. R. F. Harrington, *Field Computation by Moment Methods*, IEEE Press.
5. W. C. Gibson, *The Method of Moments in Electromagnetics*, CRC Press.
6. CST Studio Suite documentation, Microwave Studio time-domain solver and far-field monitors.
