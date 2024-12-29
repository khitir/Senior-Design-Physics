#vsim  format input:

#wz:
sqrt(w0**2 * (1 + (x**2 / (w0**2 * omega0 / (2 * c))**2)))

#betaL:
(x * gamma * delta_omega) / wz

#betaBAL:
sqrt(1 + betaL**2)

#delta_omegaL:
delta_omega / betaBAL

#omega0L:
(omega0 + betaL / betaBAL**2 * y / (w0 * delta_omega))

#dOmega0L:
omega0L - omega0

#Rz:
x * (1 + ((omega0 * w0**2) / (2 * c))**2 / x**2)

#Rinvz:
x / (x**2 + ((omega0 * w0**2) / (2 * c))**2)

#phi0L:
(Rinvz * (0**2 + (y + x * gamma * (omega0 - omega0L))**2) * omega0L) / (2 * c) + 0.5 * phi2In * (-omega0 + omega0L)**2 + (1 / c) * omega0L * (y * gamma * (-(omega0 - omega0L)) + x * (1 - 0.5 * gamma**2 * (-(omega0 - omega0L))**2))

#phi1L:
-phi2In * (omega0 - omega0L) + (1 / (2 * c)) * Rinvz * (0**2 + (y + x * gamma * (omega0 - omega0L))**2 - 2 * x * gamma * (y + x * gamma * (omega0 - omega0L)) * omega0L) + (1 / c) * ((y * gamma + x * gamma**2 * (omega0 - omega0L)) * omega0L + y * gamma * (-(omega0 - omega0L)) + x * (1 - 0.5 * gamma**2 * (-(omega0 - omega0L))**2))

#phi2L:
(1 / (2 * c)) * (2 * y * gamma - 2 * Rinvz * y * x * gamma + c * phi2In + 2 * x * gamma**2 * omega0 - 2 * Rinvz * x**2 * gamma**2 * omega0 - 3 * x * gamma**2 * omega0L + 3 * Rinvz * x**2 * gamma**2 * omega0L)

#Real_At:
(cos(omega0L * (t-b)) * w0 / wz * exp(-(y**2 / (wz**2 * betaBAL**2))) * exp(-(0**2 / wz**2)) * 1 / sqrt(1 + x**2 / ((omega0 * w0**2) / (2 * c))**2) * (exp(-((delta_omegaL**2 * ((t-b) - phi1L)**2) / (4 * (1 + delta_omegaL**4 * phi2L**2)))) * cos(phi0L - (delta_omegaL**4 * ((t-b) - phi1L)**2 * phi2L) / (4 * (1 + delta_omegaL**4 * phi2L**2)))) * 1 / (2 * sqrt(PI) * sqrt(1 / delta_omegaL**2)))

#ramp:  Vsim doesn't like this one because of & symbol
#((((t-b) > 0) && ((t-b) < tRamp)) * (sin(PI*(t-b)/(2*tRamp))**2) + ((t-b) >= tRamp)*1)


# new ramp try
((((t-b) > 0) * ((t-b) < tRamp)) * (sin(PI*(t-b)/(2*tRamp))**2) + ((t-b) >= tRamp)*1)


#Final E, is Real_at * ramp
((cos(omega0L * (t-b)) * w0 / wz * exp(-(y**2 / (wz**2 * betaBAL**2))) * exp(-(0**2 / wz**2)) * 1 / sqrt(1 + x**2 / ((omega0 * w0**2) / (2 * c))**2) * (exp(-((delta_omegaL**2 * ((t-b) - phi1L)**2) / (4 * (1 + delta_omegaL**4 * phi2L**2)))) * cos(phi0L - (delta_omegaL**4 * ((t-b) - phi1L)**2 * phi2L) / (4 * (1 + delta_omegaL**4 * phi2L**2)))) * 1 / (2 * sqrt(PI) * sqrt(1 / delta_omegaL**2))))*((((t-b) > 0) * ((t-b) < tRamp)) * (sin(PI*(t-b)/(2*tRamp))**2) + ((t-b) >= tRamp)*1)


# updated, normalized, real E function with time shift b and ramp up, vsim ready
((1 / (2 * sqrt(PI))) * ((sqrt(4*epsilon/(delta_omega * w0**2 * sqrt(PI/2)))) * w0 / wz) * cos(omega0L * (t-b) + (delta_omegaL**4 * ((t-b) - phi1L)**2 * phi2L) / (4 * (1 + delta_omegaL**4 * phi2L**2)) + atan(x/zR) - phi0L - (atan(phi2L * delta_omegaL**2)) / 2) * exp(-(y**2 / (wz**2 * betaBAL**2)) - (1/4) * (delta_omegaL**2 * ((t-b) - phi1L)**2) / (1 + (delta_omegaL**2 * phi2L)**2)) * (delta_omegaL / (1 + phi2L**2 * delta_omegaL**4)**(1/4))*((((t-b) > 0) * ((t-b) < tRamp)) * (sin(PI*(t-b)/(2*tRamp))**2) + ((t) >= tRamp)*1))
