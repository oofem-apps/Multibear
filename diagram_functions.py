import numpy as np
import math
import globvars


def ACI(Bp_y,globvar):

    fcm = globvar.concretes[globvar.fcm]
    A1 = Bp_y * globvar.Bp
    Bn = 0.85 * fcm * A1


    e = globvar.Bx / 2 - Bp_y / 2
    N_max = Bn
    M_max = Bn * e

    f_b = 0.85 * fcm * np.ones_like(Bp_y)
    d = Bp_y

    return N_max, M_max, f_b, d

def id_fcm(Bp_y,globvar):
    fcm = globvar.concretes[globvar.fcm]
    f_b = fcm * np.ones_like(Bp_y)
    d = Bp_y

    A1 = Bp_y * globvar.Bp
    e = globvar.Bx / 2 - Bp_y / 2
    N_max = A1 * f_b
    M_max = N_max * e
    return N_max, M_max, f_b, d

def ASCE(globvar):

    fcm = globvar.concretes[globvar.fcm]
    Bp_y = globvar.Bx / 2

    f_b = 2 * fcm
    A1 = Bp_y * globvar.Bp

    e = globvar.Bx / 2 - Bp_y / 2
    N_max = A1 * f_b
    M_max = N_max * e

    return N_max, M_max, f_b, Bp_y

def CTU (Bp_y,globvar):

    b_n = globvar.load_plate_width
    Bx = globvar.Bx
    d_n = Bp_y / Bx
    cover = globvar.cover


    e = Bx/2 - Bp_y/2

    diam_L = globvar.rebars_CS[globvar.DL].diam
    rebars_CS = (math.pi * diam_L ** 2) / 4
    s = globvar.s
    if globvar.flag_ignore_cover == 2:
        cover = 0
    d_l = Bx - 2 * cover

    fy = globvar.fy_vert

    ## fc
    fcm = globvar.concretes[globvar.fcm]
    fc = fcm

    ## fb_plain
    kd_plain = 4.
    kb_plain = 17.5

    Beta_b_plain = kb_plain * (1 - b_n)**2.
    Beta_d_plain = kd_plain * (1 - d_n) + d_n

    fb_plain = Beta_b_plain * Beta_d_plain

    ## fb_reinf

    # TIES
    if globvar.topology_layout == 1:

        kb_reinf = 2.9

        ke = ((1 - 0.5 * s / d_l)**2 ) / 3

        sigma_L = (2 * rebars_CS * fy) / (Bx * s)
        sigma_L_eff = sigma_L * ke
        Beta_sigma_L = (3 / sigma_L_eff)**0.3

    # MSR
    elif globvar.topology_layout == 0:

        kb_reinf = 1.4

        ke = (1 - 0.5 * s / d_l)**2

        sigma_L = (2 * rebars_CS * fy) / (d_l * s)
        sigma_L_eff = sigma_L * ke
        sigma_L_eff = 5.5
        Beta_sigma_L = 1

    # error
    else:
        print("Unknow topology")
        return 0


    Beta_b_reinf = kb_reinf * (1 - b_n)**2.

    b_n_c = 0.6
    d_n_c = 0.5

    alpha_b = np.maximum( (b_n - b_n_c) / (1 - b_n_c) , 0)
    Beta_bd_reinf = np.minimum(d_n / d_n_c , 1) * (1 - alpha_b) + d_n * alpha_b
    #fcc = Beta_sigma_L * 6.35 * sigma_L_eff**0.82 * (1 + 0.0095 * (fc - 28)**0.82) #remove staré
    fcc = Beta_sigma_L * 6.35 * sigma_L_eff**0.82 * (fc / 28)**0.2

    fb_reinf = ( fcc * (1 + Beta_b_reinf) * Beta_bd_reinf)

    ####
    #results
    ###

    plain_concrete = fc + fb_plain
    fb = fc + fb_plain + fb_reinf

    N_max = Bx**2 * fb * b_n * d_n
    M_max = N_max * e

    N_max_plain = Bx**2 * plain_concrete * b_n * d_n
    M_max_plain = N_max_plain * e

    d = Bp_y

    return N_max, M_max, fb, d, N_max_plain, M_max_plain, plain_concrete

def AS_concrete (Bp_y,globvar):

    b_n = globvar.load_plate_width
    Bx = globvar.Bx
    d_n = Bp_y / Bx


    ## fc
    fcm = globvar.concretes[globvar.fcm]
    fc = fcm

    ## fb_plain
    kd_plain = 4.
    kb_plain = 17.5

    Beta_b_plain = kb_plain * (1 - b_n)**2.
    Beta_d_plain = kd_plain * (1 - d_n) + d_n

    fb_plain = Beta_b_plain * Beta_d_plain

    plain_concrete = fc + fb_plain

    return plain_concrete


def AS_MSR (Bp_y,globvar):

    b_n = globvar.load_plate_width
    Bx = globvar.Bx
    d_n = Bp_y / Bx
    cover = globvar.cover

    e = Bx / 2 - Bp_y / 2

    diam_L = globvar.rebars_CS[globvar.DL].diam

    # if tie layout is selected
    if globvar.topology_layout == 1:
        diam_L = globvar.DL_equivalent

    rebars_CS = (math.pi * diam_L ** 2) / 4
    s = globvar.s

    #if globvar.flag_ignore_cover == 2:
    #    cover = 0

    d_l = Bx - 2 * cover

    fy = globvar.fy_vert

    ## fc
    fcm = globvar.concretes[globvar.fcm]
    fc = fcm

    ## concrete
    plain_concrete = AS_concrete(Bp_y,globvar)

    ## MSR
    kb_reinf = 1.4

    ke = (1 - 0.5 * s / d_l) ** 2

    sigma_L = (2 * rebars_CS * fy) / (d_l * s)
    sigma_L_eff = sigma_L * ke
    Beta_sigma_L = 1


    Beta_b_reinf = kb_reinf * (1 - b_n) ** 2.

    b_n_c = 0.6
    d_n_c = 0.5

    alpha_b = np.maximum((b_n - b_n_c) / (1 - b_n_c), 0)
    Beta_bd_reinf = np.minimum(d_n / d_n_c, 1) * (1 - alpha_b) + d_n * alpha_b
    fcc = Beta_sigma_L * 6.35 * sigma_L_eff ** 0.82 * (fc / 28) ** 0.2

    fb_reinf = (fcc * (1 + Beta_b_reinf) * Beta_bd_reinf)

    ####
    # results
    ###


    fb = plain_concrete + fb_reinf

    N_max = Bx ** 2 * fb * b_n * d_n
    M_max = N_max * e

    d = Bp_y

    return (N_max, M_max, fb, d)

def AS_TIE (Bp_y,globvar):

    b_n = globvar.load_plate_width
    Bx = globvar.Bx
    d_n = Bp_y / Bx
    cover = globvar.cover

    e = Bx / 2 - Bp_y / 2

    diam_L = globvar.rebars_CS[globvar.DT].diam

    # if msr layout is selected
    if globvar.topology_layout == 0:
        diam_L = globvar.DT_equivalent

    rebars_CS = (math.pi * diam_L ** 2) / 4
    s = globvar.s

    #if globvar.flag_ignore_cover == 2:
    #   cover = 0

    d_l = Bx - 2 * cover

    fy = globvar.fy_vert

    ## fc
    fcm = globvar.concretes[globvar.fcm]
    fc = fcm

    ## concrete
    plain_concrete = AS_concrete(Bp_y, globvar)

    ## TIES
    kb_reinf = 2.9

    ke = ((1 - 0.5 * s / d_l) ** 2) / 3

    sigma_L = (2 * rebars_CS * fy) / (Bx * s)
    sigma_L_eff = sigma_L * ke
    Beta_sigma_L = (3 / sigma_L_eff) ** 0.3


    Beta_b_reinf = kb_reinf * (1 - b_n) ** 2.

    b_n_c = 0.6
    d_n_c = 0.5

    alpha_b = np.maximum((b_n - b_n_c) / (1 - b_n_c), 0)
    Beta_bd_reinf = np.minimum(d_n / d_n_c, 1) * (1 - alpha_b) + d_n * alpha_b
    fcc = Beta_sigma_L * 6.35 * sigma_L_eff ** 0.82 * (fc / 28) ** 0.2

    fb_reinf = (fcc * (1 + Beta_b_reinf) * Beta_bd_reinf)

    ####
    # results
    ###


    fb = plain_concrete + fb_reinf

    N_max = Bx ** 2 * fb * b_n * d_n
    M_max = N_max * e

    d = Bp_y

    return N_max, M_max, fb, d

def AS_plain_concrete (Bp_y,globvar):

    b_n = globvar.load_plate_width
    Bx = globvar.Bx
    d_n = Bp_y / Bx

    e = Bx / 2 - Bp_y / 2

    ## concrete
    plain_concrete = AS_concrete(Bp_y, globvar)

    ####
    # results
    ###

    N_max = Bx ** 2 * plain_concrete * b_n * d_n
    M_max = N_max * e

    d = Bp_y

    return N_max , M_max , plain_concrete, d