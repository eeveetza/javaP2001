package main;



// Recommendation ITU-R P.2001

public class P2001 {
    //
    // Class implementation of Recommendation ITU-R P.2001-4
    //
    //
    //  Copyright (c) 2016-2019, Ivica Stevanovic
    //  All rights reserved.
    //
    // Redistribution and use in source and binary forms, with or without
    // modification, are permitted provided that the following conditions are
    // met:
    //
    //     * Redistributions of source code must retain the above copyright
    //       notice, this list of conditions and the following disclaimer.
    //     * Redistributions in binary form must reproduce the above copyright
    //       notice, this list of conditions and the following disclaimer in
    //       the documentation and/or other materials provided with the distribution
    //
    //
    ////
    // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    // ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    // LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    // CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    // POSSIBILITY OF SUCH DAMAGE.
    //
    // THE AUTHORS AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
    ////


    private P2001DigitalMaps maps;


    public double tl_p2001(P2001DigitalMaps maps, double[] d, double[] h, int[] z, double GHz, double Tpc, double Phire, double Phirn, double Phite,
                           double Phitn, double Hrg, double Htg, double Grx, double Gtx, int FlagVP) {
        //  tl_p2001 WRPM in the frequency range 30 MHz to 50 GHz ITU-R P.2001-3
        //   This function computes path loss due to both signal enhancements and fading
        //   over the range from 0% to 100% of an average year according to the
        //   general purpose wide-range model as described in Recommendation ITU-R
        //   P.2001-3. The model covers the frequency range from 30 MHz to 50 GHz
        //   and it is most accurate for distances from 3 km to at least 1000 km.
        //   There is no specific lower limit, although the path length
        //   must be greater than zero. A prediction of basic transmission loss less
        //   than 20 dB should be considered unreliable. Similarly, there is no
        //   specific maximum distance. Antennas heights above ground level must be
        //   greater than zero. There is no specific maximum height above ground.
        //   The method is believed to be reliable for antenna altitudes
        //   up to 8000 m above sea level.
        //
        // Inputs
        //
        // Variable    Unit    Ref     Description
        // maps                        Object containing all the Digital Maps necessary for computation
        // d           km      (2.1a)  Distance from transmitter of i-th profile point
        // h           m       (2.1b)  Height of i-th profile point (amsl)
        // z           z       (2.1c)  Zone code at distance di from transmitter
        //                             (0 = Sea, 3 = Coastal Land, 4 = Inland)
        // GHz	       GHz		T.2.2.1	Frequency
        // Tpc         %   	    T.2.2.1	Percentage of average year for which the predicted basic transmission loss is not exceeded
        // Phire       deg		T.2.2.1	Receiver longitude, positive to east
        // Phirn       deg		T.2.2.1	Receiver latitude, positive to north
        // Phite       deg		T.2.2.1	Transmitter longitude, positive to east
        // Phitn       deg		T.2.2.1	Transmitter latitude, positive to north
        // Hrg         m    	T.2.2.1	Receiving antenna height above ground
        // Htg         m    	T.2.2.1	Transmitting antenna height above ground
        // Grx         dBi		T.2.2.1	Receiving antenna gain in the direction of the ray to the transmitting antenna
        // Gtx         dBi		T.2.2.1	Transmitting antenna gain in the direction of the ray to the receiving antenna
        // FlagVp               T.2.2.1	Polarisation: 1 = vertical; 0 = horizontal
        //
        // Outputs:
        //
        // Lb         Basic transmission loss not exceeded Tpc % time.
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    19JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    29MAY17     Ivica Stevanovic, OFCOM         Corrected a bug (typo) in dl_bull_smooth
        //     v2    13JUN17     Ivica Stevanovic, OFCOM         replaced load function calls to increase computational speed
        //     v3    12JUL17     Ivica Stevanovic, OFCOM         Initial version in Java
        //     v4    19JUL17     Ivica Stevanovic, OFCOM         Introduced digital maps
        //     v5    29OCT19     Ivica Stevanovic, OFCOM         Changes in angular distance dependent loss according to ITU-R P.2001-3
        //     v6    13JUL21     Ivica Stevanovic, OFCOM         Changes in free-space loss according to ITU-R P.2001-4
        ////


        // Constants

        double c0 = 2.998e8;
        double Re = 6371;

        //// 3.1 Limited percentage time

        double Tpcp = Tpc + 0.00001 * (50 - Tpc) / 50;      // Eq (3.1.1)
        double Tpcq = 100 - Tpcp;                     // Eq (3.1.2)

        // 3.2 Path length, intermediate points, and fraction over sea

        int n = d.length;
        double dt = d[n - 1];                        // Eq (3.2.1)

        int FlagShort = 0;
        if (dt < 0.1) {
            FlagShort = 1;
        }

        // Calculate the longitude and latitude of the mid-point o the path, Phime,
        // and Phimn for dpnt = 0.5dt
        double dpnt = 0.5 * dt;
        double Phime = 0;
        double Phimn = 0;
        double Bt2r = 0;
        double Dgc = 0;


        double[] gcp = great_circle_path(Phire, Phite, Phirn, Phitn, Re, dpnt);

        Phime = gcp[0];
        Phimn = gcp[1];
        Bt2r = gcp[2];
        Dgc = gcp[3];

        // Calculate the ground height in masl at the mid-point of the profile
        // according to whether the number of profile points n is odd or even
        // (3.2.2)

        int mp = 0;
        double Hmid = 0;

        if (n % 2 == 1) { //n is odd (3.2.2a)
            mp = (n + 1) / 2 - 1;
            Hmid = h[mp];
        } else {             //n is even (3.2.2b)
            mp = n / 2 - 1;
            Hmid = 0.5 * (h[mp] + h[mp + 1]);
        }

        //Hmid = 2852.9;

        // Set the fraction of the path over sea, omega, (radio meteorological code 1)

        double omega = path_fraction(d, z, 1);

        int FlagSea = 0;
        if (omega >= 0.75) {
            FlagSea = 1;
        }

        //// 3.3 Antenna altitudes and path inclinations

        // The Tx and Rx heights masl according to (3.3.1)

        double Hts = Htg + h[0];
        double Hrs = Hrg + h[n - 1];
        double H1 = h[0];
        double Hn = h[n - 1];

        // Assign the higher and lower antenna heights above sea level (3.3.2)

        double Hhi = Math.max(Hts, Hrs);
        double Hlo = Math.min(Hts, Hrs);


        // Calculate the positive value of path inclination (3.3.3)

        double Sp = (Hhi - Hlo) / dt;

        //// 3.4 Climatic parameters
        // 3.4.1 Refractivity in the lowest 1 km


        // Find SdN from file DN_Median.txt for the path mid-pint at Phime (lon),
        // Phimn (lat) - as a bilinear interpolation

        double SdN = maps.GetSdN(Phime, Phimn);

        // Obtain Nd1km50 as in (3.4.1.1)

        double Nd1km50 = -SdN;

        // Find SdNsup from DN_SupSlope for the mid-point of the path

        double SdNsup = maps.GetSdNsup(Phime, Phimn);

        // Find SdNsub from DN_SubSlope for the mid-point of the path

        double SdNsub = maps.GetSdNsub(Phime, Phimn);

        // Obtain Nd1kmp as in (3.4.1.2)

        double Nd1kmp = 0;

        if (Tpcp < 50) {
            Nd1kmp = Nd1km50 + SdNsup * Math.log10(0.02 * Tpcp);
        } else {
            Nd1kmp = Nd1km50 - SdNsub * Math.log10(0.02 * Tpcq);
        }

        // 3.4.2 Refractivity in the lowest 65 m
        // Obtain Nd65m1 from file dndz_01.txt for the midpoint of the path

        double Nd65m1 = maps.GetNd65m1(Phime, Phimn);

        //// 3.5 Effective Earth-radius geometry
        // Median effective Earth radius (3.5.1)

        double Reff50 = 157 * Re / (157 + Nd1km50);

        // Effective Earth curvature (3.5.2)

        double Cp = (157 + Nd1kmp) / (157 * Re);

        // Effective Earth radius exceeded for p% time limited not to become
        // infinite (3.5.3)

        double Reffp = 1e6;

        if (Cp > 1e-6) {
            Reffp = 1 / Cp;
        }

        // The path length expressed as the angle subtended by d km at the center of
        // a sphere of effective Earth radius (3.5.4)

        double Thetae = dt / Reff50; // radians

        //// 3.6 Wavelength (3.6.1)

        double Wave = 1e-9 * c0 / GHz;

        //// 3.7 Path classification and terminal horizon parameters
        //  3.8 Effective heights and path roughness parameter

        double[] seh = smooth_earth_heights(d, h, Hts, Hrs, Reff50, Wave);
        double Thetat = seh[0];
        double Thetar = seh[1];
        double Thetatpos = seh[2];
        double Thetarpos = seh[3];
        double Dlt = seh[4];
        double Dlr = seh[5];
        int Ilt = (int) seh[6];
        int Ilr = (int) seh[7];
        double Hstip = seh[8];
        double Hsrip = seh[9];
        double Hstipa = seh[10];
        double Hsripa = seh[11];
        double Htea = seh[12];
        double Hrea = seh[13];
        double Mses = seh[14];
        double Hm = seh[15];
        double Hst = seh[16];
        double Hsr = seh[17];
        double Htep = seh[18];
        double Hrep = seh[19];
        int FlagLos50 = (int) seh[20];

        //// 3.9 Troposhperic-scatter path segments

        double[] trp = tropospheric_path(dt, Hts, Hrs, Thetae, Thetatpos, Thetarpos, Reff50, Phire, Phite, Phirn, Phitn, Re);

        double Dtcv = trp[0];
        double Drcv = trp[1];
        double Phicve = trp[2];
        double Phicvn = trp[3];
        double Hcv = trp[4];
        double Phitcve = trp[5];
        double Phitcvn = trp[6];
        double Phircve = trp[7];
        double Phircvn = trp[8];

        //// 3.10 Gaseous absorbtion on surface paths
        // Use the method given in Attachment F, Sec. F.2 to calculate gaseous
        // attenuations due to oxygen, and for water vapour under both non-rain and
        // rain conditions for a surface path

        double[] gas = gaseous_abs_surface(maps, Phime, Phimn, Hmid, Hts, Hrs, dt, GHz);
        double Aosur = gas[0];
        double Awsur = gas[1];
        double Awrsur = gas[2];
        double Gamo = gas[3];
        double Gamw = gas[4];
        double Gamwr = gas[5];
        double Wvsurmid = gas[6];

        double Agsur = Aosur + Awsur;                 // Eq (3.10.1)

        //// Sub-model 1

        // Calculate the diffraction loss not exceeded for p% time, as described in
        // Attachment A

        double[] dlp = dl_p(d, h, Hts, Hrs, Htep, Hrep, GHz, omega, Reffp, Cp, FlagVP);
        double Ld = dlp[0];
        double Ldsph = dlp[1];
        double Ldba = dlp[2];
        double Ldbs = dlp[3];
        double Ldbka = dlp[4];
        double Ldbks = dlp[5];
        int FlagLospa = (int) dlp[6];
        int FlagLosps = (int) dlp[7];

        //Ld = Ld_pol(FlagVP + 1);
        //Ldsph = Ldsph_pol(FlagVP + 1);

        // Use the method given in Attachemnt B.2 to calculate the notional
        // clear-air zero-fade exceedance percentage time Q0ca

        double Q0ca = multi_path_activity(GHz, dt, Hts, Hrs, Dlt, Dlr, h[Ilt], h[Ilr], Hlo, Thetat, Thetar, Sp, Nd65m1, Phimn, FlagLos50);

        // Perform the preliminary rain/wet-snow calculations in Attachment C.2
        // with the following inputs (4.1.1)

        double phi_e = Phime;
        double phi_n = Phimn;
        double h_rainlo = Hlo;
        double h_rainhi = Hhi;
        double d_rain = dt;

        double[] pfi = precipitation_fade_initial(maps, GHz, Tpcq, phi_n, phi_e, h_rainlo, h_rainhi, d_rain, FlagVP);

        double a = pfi[0];
        double b = pfi[1];
        double c = pfi[2];
        double dr = pfi[3];
        double Q0ra = pfi[4];
        double Fwvr = pfi[5];
        double kmod = pfi[6];
        double alpha_mod = pfi[7];
        int ngp = (int) pfi[8];
        double[] Gm = new double[ngp];
        double[] Pm = new double[ngp];
        for (int i = 0; i < ngp; i++) {
            Gm[i] = pfi[9 + i];
            Pm[i] = pfi[9 + ngp + i];
        }

        int flagrain = (int) pfi[9 + 2 * ngp];

        // Calculate A1 using (4.1.2)

        int flagtropo = 0; // Normal propagation close to the surface of the Earth

        double A1 = Aiter(Tpcq, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);

        // Calculate the sub-model 1 basic transmission loss not exceeded for p%  time
        double dfs = Math.sqrt(dt * dt + Math.pow((Hts - Hrs) / 1000.0, 2.0));
        //double Lbfs = tl_free_space(GHz, dt);                              // Eq (3.11.1)
        double Lbfs = tl_free_space(GHz, dfs);                              // Eq (3.11.1)

        double Lbm1 = Lbfs + Ld + A1 + Fwvr * (Awrsur - Awsur) + Agsur;      // Eq (4.1.4)

        //// Sub-model 2. Anomalous propagation

        // Use the method given in Attachment D to calculate basic transmission loss
        // not exceeded for p% time due to anomalous propagation Eq (4.2.1)

        double[] tlf = tl_anomalous_reflection(GHz, d, z, Hts, Hrs, Htea, Hrea, Hm, Thetat, Thetar, Dlt, Dlr, Phimn, omega, Reff50, Tpcp, Tpcq);

        double Lba = tlf[0];
        double Aat = tlf[1];
        double Aad = tlf[2];
        double Aac = tlf[3];
        double Dct = tlf[4];
        double Dcr = tlf[5];
        double Dtm = tlf[6];
        double Dlm = tlf[7];

        double Lbm2 = Lba + Agsur;

        //// Sub-model 3. Troposcatter propagation

        // Use the method given in Attachment E to calculate the troposcatter basic
        // transmission loss Lbs as given by equation (E.17)

        double[] tlt = tl_troposcatter(maps, GHz, dt, Thetat, Thetar, Thetae, Phicvn, Phicve, Phitn, Phite, Phirn, Phire, Gtx, Grx, Reff50, Tpcp);

        double Lbs = tlt[0];
        double Thetas = tlt[1];
        int Ztropo = (int) tlt[2];

        // To avoid under-estimating troposcatter for short paths, limit Lbs (E.17)

        Lbs = Math.max(Lbs, Lbfs);

        // Perform the preliminary rain/wet-snow calculations in Attachment C.2 from
        // the transmitter to common-volume path segment with the following inputs (4.3.1)

        phi_e = Phitcve;
        phi_n = Phitcvn;
        h_rainlo = Hts;
        h_rainhi = Hcv;
        d_rain = Dtcv;

        pfi = precipitation_fade_initial(maps, GHz, Tpcq, phi_n, phi_e, h_rainlo, h_rainhi, d_rain, FlagVP);

        a = pfi[0];
        b = pfi[1];
        c = pfi[2];
        dr = pfi[3];
        Q0ra = pfi[4];
        double Fwvrtx = pfi[5];
        kmod = pfi[6];
        alpha_mod = pfi[7];
        ngp = (int) pfi[8];
        double[] Gm1 = new double[ngp];
        double[] Pm1 = new double[ngp];
        for (int i = 0; i < ngp; i++) {
            Gm1[i] = pfi[9 + i];
            Pm1[i] = pfi[9 + ngp + i];
        }

        flagrain = (int) pfi[9 + 2 * ngp];


        // Calculate A1 using (4.1.2)

        flagtropo = 1; // for troposcatter

        double A2t = Aiter(Tpcq, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm1, Pm1, flagrain);
        //System.out.printf(     "out   =  %20.16g\n"  ,A2t);
        // Perform the preliminary rain/wet-snow calculations in Attachment C.2 from
        // the receiver to common-volume path segment with the following inputs (4.3.1)

        phi_e = Phircve;
        phi_n = Phircvn;
        h_rainlo = Hrs;
        h_rainhi = Hcv;
        d_rain = Drcv;

        pfi = precipitation_fade_initial(maps, GHz, Tpcq, phi_n, phi_e, h_rainlo, h_rainhi, d_rain, FlagVP);

        a = pfi[0];
        b = pfi[1];
        c = pfi[2];
        dr = pfi[3];
        Q0ra = pfi[4];
        double Fwvrrx = pfi[5];
        kmod = pfi[6];
        alpha_mod = pfi[7];
        ngp = (int) pfi[8];
        double[] Gm2 = new double[ngp];
        double[] Pm2 = new double[ngp];
        for (int i = 0; i < ngp; i++) {
            Gm2[i] = pfi[9 + i];
            Pm2[i] = pfi[9 + ngp + i];
        }

        flagrain = (int) pfi[9 + 2 * ngp];


        // Calculate A1 using (4.1.2)

        flagtropo = 1; // for troposcatter

        double A2r = Aiter(Tpcq, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm2, Pm2, flagrain);

        double A2 = (A2t * (1 + 0.018 * Dtcv) + A2r * (1 + 0.018 * Drcv)) / (1 + 0.018 * dt);  // Eq (4.3.6)

        // Use the method given in Attachment F.3 to calculate gaseous attenuations
        // due to oxygen and for water vapour under both non-rain and rain
        // conditions for a troposcatter path (4.3.7)

        double[] gat = gaseous_abs_tropo(maps, Phite, Phitn, Phire, Phirn, h[0], h[h.length - 1], Thetatpos, Thetarpos, Dtcv, Drcv, GHz);

        double Aos = gat[0];
        double Aws = gat[1];
        double Awrs = gat[2];
        double Aotcv = gat[3];
        double Awtcv = gat[4];
        double Awrtcv = gat[5];
        double Aorcv = gat[6];
        double Awrcv = gat[7];
        double Awrrcv = gat[8];
        double Wvsurtx = gat[9];
        double Wvsurrx = gat[10];

        // Total gaseous attenuation under non-rain conditions is given by (4.3.7)

        double Ags = Aos + Aws;

        //// Sub-model 3 basic transmission loss (4.3.8)

        double Lbm3 = Lbs + A2 + 0.5 * (Fwvrtx + Fwvrrx) * (Awrs - Aws) + Ags;

        //// 4.4 Sub-model 4. Sporadic - E

        double[] tse = tl_sporadic_e(maps, GHz, dt, Thetat, Thetar, Phimn, Phime, Phitn, Phite, Phirn, Phire, Dlt, Dlr, Reff50, Re, Tpcp);

        double Lbm4 = tse[0];
        double Lbes1 = tse[1];
        double Lbes2 = tse[2];
        double Lp1t = tse[3];
        double Lp2t = tse[4];
        double Lp1r = tse[5];
        double Lp2r = tse[6];
        double Gam1 = tse[7];
        double Gam2 = tse[8];
        double Foes1 = tse[9];
        double Foes2 = tse[10];
        double Phi1qe = tse[11];
        double Phi1qn = tse[12];
        double Phi3qe = tse[13];
        double Phi3qn = tse[14];

        //// 5 Combining sub-model results

        // 5.1 Combining sub-models 1 and 2

        double Lm = Math.min(Lbm1, Lbm2);

        double Lbm12 = Lm - 10 * Math.log10(Math.pow(10, (-0.1 * (Lbm1 - Lm))) + Math.pow(10, (-0.1 * (Lbm2 - Lm))));

        // 5.2  Combining sub-models 1+2, 3, and 4

        Lm = Math.min(Lbm12, Lbm3);
        Lm = Math.min(Lm, Lbm4);

        double Lb = Lm - 5 * Math.log10(Math.pow(10, (-0.2 * (Lbm12 - Lm))) + Math.pow(10, (-0.2 * (Lbm3 - Lm))) + Math.pow(10, (-0.2 * (Lbm4 - Lm))));
        return Lb;
        /*        double[] out = new double[127];

        out[0]	=	FlagLos50;
        out[1]	=	FlagLospa;
        out[2]	=	FlagLosps;
        out[3]  =	FlagSea;
        out[4]	=	FlagShort;
        out[5]  =	A1;
        out[6]  =	A2;
        out[7]  =	A2r;
        out[8]  =	A2t;
        out[9]  =	Aac;
        out[10] =	Aad;
        out[11] =	Aat;
        out[12] =	Ags;
        out[13] =	Agsur;
        out[14] =	Aorcv;
        out[15] =	Aos;
        out[16] =	Aosur;
        out[17] =	Aotcv;
        out[18] =	Awrcv;
        out[19] =	Awrrcv;
        out[20] =	Awrs;
        out[21] =	Awrsur;
        out[22] =	Awrtcv;
        out[23] =	Aws;
        out[24] =	Awsur;
        out[25] =	Awtcv;
        out[26] =	Bt2r;
        out[27] =	Cp;
        out[28] =	dt;
        out[29] =	Dcr;
        out[30] =	Dct;
        out[31] =	Dgc;
        out[32] =	Dlm;
        out[33] =	Dlr;
        out[34] =	Dlt;
        out[35] =	Drcv;
        out[36] =	Dtcv;
        out[37] =	Dtm;
        out[38] =	Foes1;
        out[39] =	Foes2;
        out[40] =	omega;
        out[41] =	Fwvr;
        out[42] =	Fwvrrx;
        out[43] =	Fwvrtx;
        out[44] =	Gam1;
        out[45] =	Gam2;
        out[46] =	Gamo;
        out[47] =	Gamw;
        out[48] =	Gamwr;
        out[49] =	H1;
        out[50] =	Hcv;
        out[51] =	Hhi;
        out[52] =	Hlo;
        out[53] =	Hm;
        out[54] =	Hmid;
        out[55] =	Hn;
        out[56] =	Hrea;
        out[57] =	Hrep;
        out[58] =	Hrs;
        out[59] =	Hsrip;
        out[60]=	Hsripa;
        out[61] =	Hstip;
        out[62]=	Hstipa;
        out[63] =	Htea;
        out[64] =	Htep;
        out[65] =	Hts;
        out[66] =	Lb;
        out[67] =	Lba;
        out[68] =	Lbes1;
        out[69] =	Lbes2;
        out[70] =	Lbfs;
        out[71] =	Lbm1;
        out[72] =	Lbm2;
        out[73] =	Lbm3;
        out[74] =	Lbm4;
        out[75] =	Lbs;
        out[76] =	Ld;
        out[77] =	Ldba;
        out[78] =	Ldbka;
        out[79] =	Ldbks;
        out[80] =	Ldbs;
        out[81] =	Ldsph;
        out[82] =	Lp1r;
        out[83] =	Lp1t;
        out[84] =	Lp2r;
        out[85] =	Lp2t;
        out[86] =	Mses;
        out[87] =	d.length;
        out[88] =	Nd1km50;
        out[89] =	Nd1kmp;
        out[90] =	Nd65m1;
        out[91] =	Ilr;
        out[92] =	Ilt;
        out[93] =	0;
        out[94] =	0;
        out[95] =	0;
        out[96] =	0;
        out[97] =	Phi1qe;
        out[98] =	Phi1qn;
        out[99] =	Phi3qe;
        out[100]    =	Phi3qn;
        out[101]    =	Phicve;
        out[102]    =	Phicvn;
        out[103]    =	Phime;
        out[104]    =	Phimn;
        out[105]    =	Phircve;
        out[106]    =	Phircvn;
        out[107]    =	Phitcve;
        out[108]    =	Phitcvn;
        out[109]    =	Q0ca;
        out[110]    =	Reff50;
        out[111]    =	Reffp;
        out[112]    =	Sp;
        out[113]    =	Thetae;
        out[114]    =	Thetar;
        out[115]    =	Thetarpos;
        out[116]    =	Thetas;
        out[117]    =	Thetat;
        out[118]    =	Thetatpos;
        out[119]    =	Tpcp;
        out[120]    =	Tpcq;
        out[121]    =	Math.log10(Tpcp/(100-Tpcp));
        out[122]    =	Wave;
        out[123]    =	Wvsurmid;
        out[124]    =	Wvsurrx;
        out[125]    =	Wvsurtx;
        out[126]    =	Ztropo;

        return out; */
    }

    public double path_fraction(double[] d, int[] zone, int zone_r) {
        //path_fraction Path fraction belonging to a given zone_r
        //     omega = path_fraction(d, zone, zone_r)
        //     This function computes the path fraction belonging to a given zone_r
        //     of the great-circle path (km)
        //
        //     Input arguments:
        //     d       -   vector of distances in the path profile
        //     zone    -   vector of zones in the path profile
        //     zone_r  -   reference zone for which the fraction is computed
        //
        //     Output arguments:
        //     omega   -   path fraction belonging to the given zone_r
        //
        //     Example:
        //     omega = path_fraction(d, zone, zone_r)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    17NOV16     Ivica Stevanovic, OFCOM         First implementation in Java

        double dm = 0;

        int n = d.length;
        double delta;

        for (int i = 0; i < n; i++) {
            if (zone[i] == zone_r) {
                if (i == 0) {
                    delta = (d[1] - d[0]) / 2.0;
                } else if (i == (n - 1)) {
                    delta = (d[n - 1] - d[n - 2]) / 2.0;
                } else {
                    delta = (d[i + 1] - d[i - 1]) / 2.0;
                }

                dm = dm + delta;
            }
        }


        return dm / (d[n - 1] - d[0]);
    }


    public double[] great_circle_path(double Phire, double Phite, double Phirn, double Phitn, double Re, double dpnt) {
        //great_circle_path Great-circle path calculations according to Attachment H
        //   This function computes the great-circle intermediate points on the
        //   radio path as defined in ITU-R P.2001-3 Attachment H
        //
        //     Input parameters:
        //     Phire   -   Receiver longitude, positive to east (deg)
        //     Phite   -   Transmitter longitude, positive to east (deg)
        //     Phirn   -   Receiver latitude, positive to north (deg)
        //     Phitn   -   Transmitter latitude, positive to north (deg)
        //     Re      -   Average Earth radius (km)
        //     dpnt    -   Distance from the transmitter to the intermediate point (km)
        //
        //     Output parameters:
        //     Phipnte -   Longitude of the intermediate point (deg)
        //     Phipntn -   Latitude of the intermediate point (deg)
        //     Bt2r    -   Bearing of the great-circle path from Tx towards the Rx (deg)
        //     dgc     -   Great-circle path length (km)
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    12JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUL17     Ivica Stevanovic, OFCOM         Initial version in Java


        //// H.2 Path length and bearing

        // Difference (deg) in longitude between the terminals (H.2.1)

        double Dlon = Phire - Phite;

        // Calculate quantity r (H.2.2)

        double r = sind(Phitn) * sind(Phirn) + cosd(Phitn) * cosd(Phirn) * cosd(Dlon);

        // Calculate the path length as the angle subtended at the center of
        // average-radius Earth (H.2.3)

        double Phid = Math.acos(r);  // radians

        // Calculate the great-circle path length (H.2.4)

        double dgc = Phid * Re;  // km

        // Calculate the quantity x1 (H.2.5a)

        double x1 = sind(Phirn) - r * sind(Phitn);

        // Calculate the quantity y1 (H.2.5b)

        double y1 = cosd(Phitn) * cosd(Phirn) * sind(Dlon);

        // Calculate the bearing of the great-circle path for Tx to Rx (H.2.6)

        double Bt2r = Phire;

        if (Math.abs(x1) < 1e-9 && Math.abs(y1) < 1e-9) {
            Bt2r = Phire;
        } else {
            Bt2r = atan2d(y1, x1);
        }

        //// H.3 Calculation of intermediate path point

        // Calculate the distance to the point as the angle subtended at the center
        // of average-radius Earth (H.3.1)

        double Phipnt = dpnt / Re;  //radians

        // Calculate quantity s (H.3.2)

        double s = sind(Phitn) * Math.cos(Phipnt) + cosd(Phitn) * Math.sin(Phipnt) * cosd(Bt2r);

        // The latitude of the intermediate point is now given by (H.3.3)

        double Phipntn = asind(s); // degs

        // Calculate the quantity x2 (H.3.4a)

        double x2 = Math.cos(Phipnt) - s * sind(Phitn);

        // Calculate the quantity y2 (H.3.4b)

        double y2 = cosd(Phitn) * Math.sin(Phipnt) * sind(Bt2r);

        // Calculate the longitude of the intermediate point Phipnte (H.3.5)

        double Phipnte = Bt2r;

        if (x2 < 1e-9 && y2 < 1e-9) {
            Phipnte = Bt2r;
        } else {
            Phipnte = Phite + atan2d(y2, x2);
        }

        double[] out = new double[4];

        out[0] = Phipnte;
        out[1] = Phipntn;
        out[2] = Bt2r;
        out[3] = dgc;

        return out;
    }

    public double sind(double theta_deg) {

        double theta_rad = theta_deg * Math.PI / 180.0;
        double y = Math.sin(theta_rad);
        return y;
    }

    public double cosd(double theta_deg) {

        double theta_rad = theta_deg * Math.PI / 180.0;
        double y = Math.cos(theta_rad);
        return y;
    }

    public double atan2d(double y, double x) {

        double res = Math.atan2(y, x) * 180.0 / Math.PI;
        return res;

    }

    public double asind(double y) {

        double x = Math.asin(y) * 180.0 / Math.PI;
        return x;

    }

    public double[] smooth_earth_heights(double[] d, double[] h, double hts, double hrs, double ae, double lam) {
        //smooth_earth_heights smooth-Earth effective antenna heights according to ITU-R P.2001-3
        // [theta_t, theta_r, theta_tpos, theta_rpos, dlt, dlr, lt, lr, hstip, hsrip, hstipa, hsripa, htea, hrea, hm, hst, hsr, htep, hrep] = smooth_earth_heights(d, h, hts, hrs, ae, lam)
        // This function derives smooth-Earth effective antenna heights according to
        // Sections 3.7 and 3.8 of Recommendation ITU-R P.2001-3
        //
        // Input parameters:
        // d         -   vector of terrain profile distances from Tx [0,dtot] (km)
        // h         -   vector of terrain profile heights amsl (m)
        // hts, hrs  -   Tx and Rx antenna heights above means sea level (m)
        //               hts = htg + h(1), hrs = hrg + h(end)
        // ae        -   median effective Earth's radius
        // lam       -   wavelength (m)
        //
        // Output parameters:
        //
        // theta_t      -   Interfering antenna horizon elevation angle (mrad)
        // theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
        // theta_tpos   -   Interfering antenna horizon elevation angle limited to be positive (mrad)
        // theta_rpos   -   Interfered-with antenna horizon elevation angle limited to be positive (mrad)
        // dlt          -   Tx antenna horizon distance (km)
        // dlr          -   Rx antenna horizon distance (km)
        // lt           -   Index i in the path profile for which dlt=d(lt)
        // lr           -   Index i in tha path profile for which dlr = d(lr)
        // hstip, hrip  -   Initial smooth-surface height at Tx/Rx
        // hstipa, hripa-   Smooth surface height at Tx/Rx not exceeding ground level
        // htea, htea   -   Effective Tx and Rx antenna heigts above the smooth-Earth surface amsl for anomalous propagation (m)
        // mses         -   Smooth surface slope (m/km)
        // hm           -   The terrain roughness parameter (m)
        // hst, hsr     -   Heights of the smooth surface at the Tx and Rx ends of the paths
        // htep, hrep   -   Effective Tx and Rx antenna heights for the
        //                  spherical-earth and the smooth-pofile version of the
        //                  Bullingtong diffraction model (m)
        // FlagLos50    -   1 = Line-of-sight 50// time, 0 = otherwise
        //
        //
        // Rev   Date        Author                          Description
        // -------------------------------------------------------------------------------
        // v0    15JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab (P.452)
        // v1    15JUN16     Ivica Stevanovic, OFCOM         Modifications related to LoS path (P.452)
        // v3    15JUN16     Ivica Stevanovic, OFCOM         Initial version for P.1812
        // v4    13JUL16     Ivica Stevanovic, OFCOM         Initial version for P.2001
        // v5    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        int n = d.length;

        double dtot = d[n - 1];

        //// 3.7 Path classification and terminal horizon parameters

        // Highest elevation angle to an intermediate profile point, relative to the
        // horizontal at the transmitter (3.7.1)

        double theta_tim = ((h[1] - hts) / d[1] - 500 * d[1] / ae);

        for (int ii = 2; ii < n - 1; ii++) {
            theta_tim = Math.max(theta_tim, ((h[ii] - hts) / d[ii] - 500 * d[ii] / ae));
        }


        // Elevation angle of the receiver as viewed by the transmitter, assuming a
        // LoS path (3.7.2)

        double theta_tr = (hrs - hts) / dtot - 500 * dtot / ae;
        int FlagLos50 = 1;
        int kindex = 1;
        int lt = 0;
        double dlt = 0;
        double dlr = 0;
        int lr = 0;
        double theta_t = 0;
        double theta_r = 0;

        if (theta_tim < theta_tr) { // path is LoS
            FlagLos50 = 1;
            int ii = 1;
            double numax = (h[ii] + 500 * d[ii] * (dtot - d[ii]) / ae - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot) *
                    Math.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])));             // Eq (3.7.3)

            for (ii = 2; ii < n - 1; ii++) {
                double nu = (h[ii] + 500 * d[ii] * (dtot - d[ii]) / ae - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot) *
                        Math.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])));             // Eq (3.7.3)
                if (nu > numax) {
                    numax = nu;
                    kindex = ii;

                }
            }

            lt = kindex;     //in order to map back to path d indices, as theta takes path indices 2 to n-1,
            dlt = d[lt];                                // Eq (3.7.4a)
            dlr = dtot - dlt;                           // Eq (3.7.4b)
            lr = lt;                                      // Eq (3.7.4d)

            theta_t = theta_tr;                        // Eq (3.7.5a)
            theta_r = -theta_tr - 1000 * dtot / ae;         // Eq (3.7.5b)

        } else {
            FlagLos50 = 0;

            // Transmitter horizon distance and profile index of the horizon point (3.7.6)
            for (int ii = 1; ii < n - 1; ii++) {
                double theta_ti = ((h[ii] - hts) / d[ii] - 500 * d[ii] / ae);
                if (theta_ti == theta_tim) {
                    lt = ii;     //in order to map back to path d indices, as theta takes path indices 2 to n-1,
                    dlt = d[lt];     // Eq (3.7.6a)
                    break;
                }
            }

            // Transmitter horizon elevation angle reltive to its local horizontal (3.7.7)
            theta_t = theta_tim;

            // Find the heighest elevation angle to an intermediate profile point,
            // relative to the horizontal at the receiver (3.7.8)
            int ii = 1;
            double theta_ri = ((h[ii] - hrs) / (dtot - d[ii]) - 500 * (dtot - d[ii]) / ae);
            double theta_rim = theta_ri;

            for (ii = 2; ii < n - 1; ii++) {
                theta_ri = ((h[ii] - hrs) / (dtot - d[ii]) - 500 * (dtot - d[ii]) / ae);
                if (theta_ri > theta_rim) {
                    theta_rim = theta_ri;
                    kindex = ii;
                }
            }

            lr = kindex;     //in order to map back to path d indices, as theta takes path indices 2 to n-1,
            dlr = dtot - d[lr];                           // Eq (3.7.9)

            // receiver horizon elevation angle relative to its local horizontal
            theta_r = theta_rim;                        // Eq (3.7.10)
        }

        // Calculate the horizon elevation angles limited such that they are
        // positive

        double theta_tpos = Math.max(theta_t, 0);                   // Eq (3.7.11a)
        double theta_rpos = Math.max(theta_r, 0);                   // Eq (3.7.11b)

        //// 3.8 Effective heights and path roughness parameter

        double v1 = 0;
        for (int ii = 1; ii < n; ii++) {
            v1 = v1 + (d[ii] - d[ii - 1]) * (h[ii] + h[ii - 1]);  // Eq (85)
        }
        double v2 = 0;

        for (int ii = 1; ii < n; ii++) {
            v2 = v2 + (d[ii] - d[ii - 1]) * (h[ii] * (2 * d[ii] + d[ii - 1]) + h[ii - 1] * (d[ii] + 2 * d[ii - 1]));  // Eq (86)
        }


        double hstip = (2 * v1 * dtot - v2) / (dtot * dtot);       // Eq (3.8.3a)
        double hsrip = (v2 - v1 * dtot) / (dtot * dtot);          // Eq (3.8.3b)

        // Smooth-surface heights limited not to exceed ground level at either Tx or Rx

        double hstipa = Math.min(hstip, h[0]);               // Eq (3.8.4a)
        double hsripa = Math.min(hsrip, h[n - 1]);             // Eq (3.8.4b)

        // The slope of the least-squares regression fit (3.8.5)

        double mses = (hsripa - hstipa) / dtot;

        // effective heights of Tx and Rx antennas above the smooth surface (3.8.6)

        double htea = hts - hstipa;
        double hrea = hrs - hsripa;

        // Path roughness parameter (3.8.7)

        int ii = lt;
        double hm = h[ii] - (hstipa + mses * d[ii]);

        for (ii = lt + 1; ii <= lr; ii++) {

            hm = Math.max(hm, h[ii] - (hstipa + mses * d[ii]));

        }

        // Smooth-surface heights for the diffraction model

        ii = 1;

        double hobs = h[ii] - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot;  //  Eq (3.8.8d)
        double alpha_obt = hobs / d[ii]; // Eq (3.8.8b)

        double alpha_obr = hobs / (dtot - d[ii]); // Eq (3.8.8c)
        for (ii = 2; ii < n - 1; ii++) {
            double HH = h[ii] - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot;

            hobs = Math.max(hobs, HH);                 // Eq (3.8.8a)

            alpha_obt = Math.max(alpha_obt, HH / d[ii]); // Eq (3.8.8b)

            alpha_obr = Math.max(alpha_obr, HH / (dtot - d[ii])); // Eq (3.8.8c)
        }
        // Calculate provisional values for the Tx and Rx smooth surface heights

        double gt = alpha_obt / (alpha_obt + alpha_obr);         // Eq (3.8.9e)
        double gr = alpha_obr / (alpha_obt + alpha_obr);         // Eq (3.8.9f)

        double hst = 0;
        double hsr = 0;

        if (hobs <= 0) {
            hst = hstip;                                // Eq (3.8.9a)
            hsr = hsrip;                                // Eq (3.8.9b)
        } else {
            hst = hstip - hobs * gt;                      // Eq (3.8.9c)
            hsr = hsrip - hobs * gr;                      // Eq (3.8.9d)
        }

        // calculate the final values as required by the diffraction model

        if (hst >= h[0]) {
            hst = h[0];                                // Eq (3.8.10a)
        }

        if (hsr > h[n - 1]) {
            hsr = h[n - 1];                              // Eq (3.8.10b)
        }

        // The terminal effective heigts for the ducting/layer-reflection model

        double htep = hts - hst;                            // Eq (3.8.11a)
        double hrep = hrs - hsr;                            // Eq (3.8.11b)

        double[] seh = new double[21];//declaration and instantiation.
        seh[0] = theta_t;
        seh[1] = theta_r;
        seh[2] = theta_tpos;
        seh[3] = theta_rpos;
        seh[4] = dlt;
        seh[5] = dlr;
        seh[6] = lt;
        seh[7] = lr;
        seh[8] = hstip;
        seh[9] = hsrip;
        seh[10] = hstipa;
        seh[11] = hsripa;
        seh[12] = htea;
        seh[13] = hrea;
        seh[14] = mses;
        seh[15] = hm;
        seh[16] = hst;
        seh[17] = hsr;
        seh[18] = htep;
        seh[19] = hrep;
        seh[20] = (double) FlagLos50;

        return seh;
    }

    public double[] tropospheric_path(double dt, double hts, double hrs, double theta_e, double theta_tpos, double theta_rpos, double ae, double phi_re, double phi_te, double phi_rn, double phi_tn, double Re) {
        //trophospheric path segments according to ITU-R P.2001-3
        // This function computes tropospheric path segments as described in Section
        // 3.9 of Recommendation ITU-R P.2001-3
        //
        // Input parameters:
        // dt        -   Path length (km)
        // hts, hrs  -   Tx/Rx antenna heights above means sea level (m)
        // theta_e   -   Angle subtended by d km at the center of a sphere of effective earth radius (rad)
        // theta_tpos-   Interfering antenna horizon elevation angle limited to be positive (mrad)
        // theta_rpos-   Interfered-with antenna horizon elevation angle limited to be positive (mrad)
        //               hts = htg + h(1)
        // ae        -   median effective Earth's radius
        // phi_re    -   Receiver longitude, positive to east (deg)
        // phi_te    -   Transmitter longitude, positive to east (deg)
        // phi_rn    -   Receiver latitude, positive to north (deg)
        // phi_tn    -   Transmitter latitude, positive to north (deg)
        // Re        -   Average Earth radius (km)
        //
        // Output parameters:
        // d_tcv     -   Horizontal path length from transmitter to common volume (km)
        // d_rcv     -   Horizontal path length from common volume to receiver (km)
        // phi_cve   -   Longitude of the common volume
        // phi_cvn   -   Latitude of the common volume
        // h_cv      -   Height of the troposcatter common volume (masl)
        // phi_tcve  -   Longitude of midpoint of the path segment from Tx to common volume
        // phi_tcvn  -   Latitude of midpoint of the path segment from Tx to common volume
        // phi_rcve  -   Longitude of midpoint of the path segment from common volume to Rx
        // phi_rcvn  -   Latitude of midpoint of the path segment from common volumen to Rx
        //
        //
        // Rev   Date        Author                          Description
        // -------------------------------------------------------------------------------
        // v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        // v1    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // Horizontal path lenght from transmitter to common volumne (3.9.1a)

        double d_tcv = (dt * Math.tan(0.001 * theta_rpos + 0.5 * theta_e) - 0.001 * (hts - hrs)) /
                (Math.tan(0.001 * theta_tpos + 0.5 * theta_e) + Math.tan(0.001 * theta_rpos + 0.5 * theta_e));

        // Limit d_tcv such that 0 <= dtcv <= dt

        if (d_tcv < 0) {
            d_tcv = 0;
        }
        if (d_tcv > dt) {
            d_tcv = dt;
        }

        // Horizontal path length from common volume to receiver (3.9.1b)

        double d_rcv = dt - d_tcv;

        // Calculate the longitude and latitude of the common volumne from the
        // transmitter and receiver longitudes and latitudes using the great circle
        // path method of Attachment H by seting d_pnt = d_tcv

        double[] gcp = great_circle_path(phi_re, phi_te, phi_rn, phi_tn, Re, d_tcv);
        double phi_cve = gcp[0];
        double phi_cvn = gcp[1];


        // Calculate the height of the troposcatter common volume (3.9.2)

        double h_cv = hts + 1000 * d_tcv * Math.tan(0.001 * theta_tpos) + 1000 * d_tcv * d_tcv / (2 * ae);

        // Calculate the longitude and latitude of the midpoint of hte path segment
        // from transmitter to common volume by setting dpnt = 0.5dtcv

        double d_pnt = 0.5 * d_tcv;

        gcp = great_circle_path(phi_re, phi_te, phi_rn, phi_tn, Re, d_pnt);
        double phi_tcve = gcp[0];
        double phi_tcvn = gcp[1];

        // Calculate the longitude and latitude of the midpoint of the path segment
        // from receiver to common volume by setting dpnt = dt - 0.5drcv

        d_pnt = dt - 0.5 * d_rcv;

        gcp = great_circle_path(phi_re, phi_te, phi_rn, phi_tn, Re, d_pnt);
        double phi_rcve = gcp[0];
        double phi_rcvn = gcp[1];


        double[] out = new double[9];
        out[0] = d_tcv;
        out[1] = d_rcv;
        out[2] = phi_cve;
        out[3] = phi_cvn;
        out[4] = h_cv;
        out[5] = phi_tcve;
        out[6] = phi_tcvn;
        out[7] = phi_rcve;
        out[8] = phi_rcvn;

        return out;

    }


    public double[] gaseous_abs_surface(P2001DigitalMaps maps, double phi_me, double phi_mn, double h_mid, double hts, double hrs, double dt, double f) {
        //specific_sea_level_attenuation Specific sea-level attenuations (Attachment F.6)
        //   This function computes specific sea-level attenuations due to oxigen
        //   and water-vapour as defined in ITU-R P.2001-3 Attachment F.6
        //   The formulas are valid for frequencies not greater than 54 GHz.
        //
        //   Input parameters:
        //   phi_me    -   Longitude of the mid-point of the path (deg)
        //   phi_mn    -   Latitude of the mid-point of the path (deg)
        //   h_mid     -   Ground height at the mid-point of the profile (masl)
        //   hts, hrs  -   Tx and Rx antenna heights above means sea level (m)
        //                 hts = htg + h(1), hrs = hrg + h(end)
        //   f         -   Frequency (GHz), not greater than 54 GHz
        //
        //   Output parameters:
        //   Aosur     -   Attenuation due to oxygen (dB)
        //   Awsur     -   Attenuation due to water-vapour under non-rain conditions (dB)
        //   Awrsur    -   Attenuation due to water-vapour under rain conditions (dB)
        //   gamma_o  -   Specific attenuation due to oxygen (dB/km)
        //   gamma_w  -   Specific attenuation due to water-vapour non-rain conditions (dB/km)
        //   gamma_wr -   Specific attenuation due to water-vapour rain conditions (dB/km)
        //   rho_sur  -   Surface water-vapour content (g/m^3)
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUN17     Ivica Stevanovic, OFCOM         replaced load function calls to increase computational speed
        //     v2    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // Obtain surface water-vapour density under non-rain conditions at the
        // midpoint of the path from the data file surfwv_50_fixed.txt

        //surfwv_50_fixed = load('DigitalMaps/surfwv_50_fixed.txt');
        double rho_sur = maps.GetRhoSur(phi_me, phi_mn);

        double h_sur = h_mid;

        // Use equation (F.6.2) to calculate the sea-level specific attenuation due
        // to water vapour under non-rain conditions gamma_w

        double[] gamma = specific_sea_level_attenuation(f, rho_sur, h_sur);
        double gamma_o = gamma[0];
        double gamma_w = gamma[1];

        // Use equation (F.5.1) to calculate the surface water-vapour density under
        // rain conditions rho_surr

        double rho_surr = water_vapour_density_rain(rho_sur, h_sur);

        // Use equation (F.6.2) to calculate the sea-level specifica ttenuation due
        // to water vapour undr rain conditions gamma_wr

        double[] gamma_xx = specific_sea_level_attenuation(f, rho_surr, h_sur);
        double gamma_wr = gamma_xx[1];


        // Calculate the height for water-vapour density (F.2.1)

        double h_rho = 0.5 * (hts + hrs);

        // Attenuation due to oxygen (F.2.2a)

        double Aosur = gamma_o * dt * Math.exp(-h_rho / 5000);

        // Attenuation due to water-vapour under non-rain conditions (F.2.2b)

        double Awsur = gamma_w * dt * Math.exp(-h_rho / 2000);

        // Attenuation due to water-vapour under non-rain conditions (F.2.2b)

        double Awrsur = gamma_wr * dt * Math.exp(-h_rho / 2000);

        double[] out = new double[7];

        out[0] = Aosur;
        out[1] = Awsur;
        out[2] = Awrsur;
        out[3] = gamma_o;
        out[4] = gamma_w;
        out[5] = gamma_wr;
        out[6] = rho_sur;

        return out;
    }


    public double[] specific_sea_level_attenuation(double f, double rho_sur, double h_sur) {
        //specific_sea_level_attenuation Specific sea-level attenuations (Attachment F.6)
        //   This function computes specific sea-level attenuations due to oxigen
        //   and water-vapour as defined in ITU-R P.2001-3 Attachment F.6
        //   The formulas are valid for frequencies not greater than 54 GHz.
        //
        //     Input parameters:
        //     f       -   Frequency (GHz), not greater than 54 GHz
        //     rho_sur -   Surface water-vapour density under non-rain conditions (g/m^3)
        //     h_sur   -   Terrain height (masl)
        //
        //     Output parameters:
        //     gamma_o -   Sea-level specific attenuation due to oxigen (dB/km)
        //     gamma_w -   Sea-level specific attenuation due to water vapour (dB/km)

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    12JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        double rho_sea = rho_sur * Math.exp(h_sur / 2000);      // Eq (F.6.2b)

        double eta = 0.955 + 0.006 * rho_sea;             // Eq (F.6.2a)

        double gamma_o = (7.2 / (f * f + 0.34) + 0.62 / (Math.pow((54 - f), 1.16) + 0.83)) * f * f * 1e-3;    // Eq (F.6.1)

        double gamma_w = (0.046 + 0.0019 * rho_sea + 3.98 * eta / (Math.pow((f - 22.235), 2) + 9.42 * eta * eta) * (1 + Math.pow(((f - 22) / (f + 22)), 2))) * f * f * rho_sea * 1e-4;    //Eq (F.6.2)

        double[] out = new double[2];
        out[0] = gamma_o;
        out[1] = gamma_w;

        return out;
    }

    public double water_vapour_density_rain(double rho_sur, double h_sur) {
        //water_vapour_density_rain Atmospheric water-vapour density in rain (Attachment F.5)
        //   This function computes atmosphoric water-vapour density in rain
        //   as defined in ITU-R P.2001-3 Attachment F.5
        //
        //     Input parameters:
        //     rho_sur -   Surface water-vapour density under non-rain conditions (g/m^3)
        //     h_sur   -   Terrain height (masl)
        //
        //     Output parameters:
        //     rho_surr-   Atmospheric water-vapour density in rain

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        double rho_surr = 0;

        if (h_sur <= 2600) {
            rho_surr = rho_sur + 0.4 + 0.0003 * h_sur;
        } else {
            rho_surr = rho_sur + 5 * Math.exp(-h_sur / 1800);
        }

        return rho_surr;
    }


    public double[] dl_p(double[] d, double[] h, double hts, double hrs, double hte, double hre, double f, double omega, double ap, double Cp, int polHV) {
        //dl_p Diffraction loss model not exceeded for p// of time according to P.2001-3
        //   [Ld, Ldsph, Ldba, Ldbs, Ldbka, Ldbks] = dl_p( d, h, hts, hrs, hte, hre, f, omega, ap, Cp )
        //
        //   This function computes the diffraction loss not exceeded for p// of time
        //   as defined in ITU-R P.2001-3 (Attachment A)
        //
        //     Input parameters:
        //     d       -   vector of distances di of the i-th profile point (km)
        //     h       -   vector hi of heights of the i-th profile point (meters
        //                 above mean sea level).
        //                 Both vectors h and d contain n+1 profile points
        //     hts     -   transmitter antenna height in meters above sea level (i=0)
        //     hrs     -   receiver antenna height in meters above sea level (i=n)
        //     hte     -   Effective height of interfering antenna (m amsl)
        //     hre     -   Effective height of interfered-with antenna (m amsl)
        //     f       -   frequency expressed in GHz
        //     omega   -   the fraction of the path over sea
        //     ap      -   Effective Earth radius (km)
        //     Cp      -   Effective Earth curvature
        //     polHV   -   polarization (0) horizontal (1) vertical
        //     Output parameters:
        //     Ldp    -   diffraction loss for the general path not exceeded for p // of the time
        //                according to Attachment A of ITU-R P.2001-3.
        //                Ldp(1) is for the horizontal polarization
        //                Ldp(2) is for the vertical polarization
        //     Ldshp  -   Spherical-Earth diffraction loss diffraction (A.2) for the actual path d and modified antenna heights
        //     Lba    -   Bullington diffraction loss for the actual path profile as calculated in A.4
        //     Lbs    -   Bullingtong diffraction loss for a smooth path profile as calculated in A.5
        //     Ldbka  -   Knife-edge diffraction loss for Bullington point: actual path
        //     Ldbks  -   Knife-edge diffraction loss for Bullington point: smooth path
        //     FlagLospa - 1 = LoS p// time for actual path, 0 = otherwise
        //     FlagLosps - 1 = LoS p// time for smooth path, 0 = otherwise

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version (P.452)
        //     v1    06JUL16     Ivica Stevanovic, OFCOM         Modifications according to P.1812
        //     v2    13JUL16     Ivica Stevanovic, OFCOM         Modifications according to P.2001
        //     v3    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version


        ////
        int n = d.length;
        double dtot = d[n - 1];

        double[] Ldsph = dl_se(dtot, hte, hre, ap, f, omega);

        double[] dba = dl_bull_actual(d, h, hts, hrs, Cp, f);
        double Ldba = dba[0];
        double Ldbka = dba[1];
        double FlagLospa = dba[2];

        double[] dbs = dl_bull_smooth(d, h, hte, hre, ap, f);
        double Ldbs = dbs[0];
        double Ldbks = dbs[1];
        double FlagLosps = dbs[2];

        double Ld = Ldba + Math.max(Ldsph[polHV] - Ldbs, 0);                  // Eq (A.1.1)

        double[] out = new double[8];
        out[0] = Ld;
        out[1] = Ldsph[polHV];
        out[2] = Ldba;
        out[3] = Ldbs;
        out[4] = Ldbka;
        out[5] = Ldbks;
        out[6] = FlagLospa;
        out[7] = FlagLosps;
        return out;
    }


    public double[] dl_se(double d, double hte, double hre, double ap, double f, double omega) {
        //dl_se spherical-Earth diffraction loss exceeded for p// time according to ITU-R P.2001-3
        //   This function computes the Spherical-Earth diffraction loss not exceeded
        //   for p// time for antenna heights hte and hre (m)
        //   as defined in Attachment A.2 of ITU-R P.2001-3
        //
        //     Input parameters:
        //     d       -   Great-circle path distance (km)
        //     hte     -   Effective height of interfering antenna (m)
        //     hre     -   Effective height of interfered-with antenna (m)
        //     ap      -   the effective Earth radius in kilometers
        //     f       -   Frequency (GHz)
        //     omega   -   the fraction of the path over sea
        //
        //     Output parameters:
        //     Ldsph   -   The spherical-Earth diffraction loss not exceeded for p// time
        //                 Ldsph(1) is for the horizontal polarization
        //                 Ldsph(2) is for the vertical polarization
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    23DEC15     Ivica Stevanovic, OFCOM         Initial version (P.452)
        //     v1    01FEB16     Ivica Stevanovic, OFCOM         Introduced dl_se_ft (P. 452)
        //     v2    06JUL16     Ivica Stevanovic, OFCOM         Initial version (P.1812)
        //     v2    13JUL16     Ivica Stevanovic, OFCOM         Initial version (P.2001)
        //     v3    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // Initialize
        double[] Ldsph = new double[]{0, 0};

        // Wavelength
        double c0 = 2.998e8;
        double lam = 1e-9 * c0 / f;

        // Calculate the marginal LoS distance for a smooth path

        double dlos = Math.sqrt(2 * ap) * (Math.sqrt(0.001 * hte) + Math.sqrt(0.001 * hre));    // Eq (A.2.1)

        if (d >= dlos) {
            // calculate diffraction loss Ldft using the method in Sec. A.3 for
            // adft = ap and set Ldsph to Ldft

            Ldsph = dl_se_ft(d, hte, hre, ap, f, omega);
            return Ldsph;
        } else {
            // calculate the smallest clearance between the curved-Earth path and
            // the ray between the antennas, hse

            double c = (hte - hre) / (hte + hre);        // Eq (A.2.2d)
            double m = 250 * d * d / (ap * (hte + hre));        // Eq (A.2.2e)

            double b = 2.0 * Math.sqrt((m + 1.0) / (3.0 * m)) * Math.cos(Math.PI / 3.0 + 1.0 / 3.0 * Math.acos(3.0 * c / 2.0 * Math.sqrt(3.0 * m / Math.pow((m + 1.0), 3.0))));   // Eq (A.2.2c)

            double dse1 = d / 2.0 * (1 + b);           // Eq (A.2.2a)
            double dse2 = d - dse1;            // Eq (A.2.2b)

            double hse = (hte - 500.0 * dse1 * dse1 / ap) * dse2 + (hre - 500.0 * dse2 * dse2 / ap) * dse1;
            hse = hse / d;                // Eq (A.2.2)

            // Calculate the required clearance for zero diffraction loss

            double hreq = 17.456 * Math.sqrt(dse1 * dse2 * lam / d);     // Eq (A.2.3)

            if (hse > hreq) {
                Ldsph[0] = 0;
                Ldsph[1] = 0;
                return Ldsph;
            } else {

                // calculate the modified effective Earth radius aem, which gives
                // marginal LoS at distance d

                double aem = 500 * Math.pow((d / (Math.sqrt(hte) + Math.sqrt(hre))), 2);     // Eq (A.2.4)

                // Use the method in Sec. A3 for adft = aem to obtain Ldft

                double[] Ldft = dl_se_ft(d, hte, hre, aem, f, omega);

                if (Ldft[0] < 0) {
                    Ldsph[0] = 0;
                    Ldsph[1] = 0;
                    return Ldsph;
                } else {
                    Ldsph[0] = (1 - hse / hreq) * Ldft[0];     // Eq (A.2.5)
                    Ldsph[1] = (1 - hse / hreq) * Ldft[1];     // Eq (A.2.5)

                }
            }
        }
        //System.out.printf(     "Ldft[0]    =  %g\n"  , Ldsph[0]);
        //System.out.printf(     "Ldft[0]    =  %g\n"  , Ldsph[1]);
        return Ldsph;
    }

    public double[] dl_se_ft(double d, double hte, double hre, double adft, double f, double omega) {
        //dl_se_ft First-term part of spherical-Earth diffraction according to ITU-R P.2001-3
        //   This function computes the first-term part of Spherical-Earth diffraction
        //   as defined in Sec. A.3 of the ITU-R P.2001-3
        //
        //     Input parameters:
        //     d       -   Great-circle path distance (km)
        //     hte     -   Effective height of interfering antenna (m)
        //     hre     -   Effective height of interfered-with antenna (m)
        //     adft    -   effective Earth radius (km)
        //     f       -   Frequency (GHz)
        //     omega   -   fraction of the path over sea
        //
        //     Output parameters:
        //     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p// time
        //                Ldft(1) is for the horizontal polarization
        //                Ldft(2) is for the vertical polarization
        //
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation (P.452)
        //     v1    06JUL16     Ivica Stevanovic, OFCOM         First implementation (P.1812)
        //     v2    13JUL16     Ivica Stevanovic, OFCOM         First implementation (P.2001)
        //     v3    05AUG17     Ivica Stevanovic, OFCOM         epsr and sigma were swapped


        ////

        // First-term part of the spherical-Earth diffraction loss over land

        double epsr = 22;
        double sigma = 0.003;

        double[] Ldft_land = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f);

        // First-term part of the spherical-Earth diffraction loss over sea

        epsr = 80;
        sigma = 5;

        double[] Ldft_sea = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f);


        // First-term spherical diffraction loss

        double[] Ldft = new double[2];
        Ldft[0] = omega * Ldft_sea[0] + (1 - omega) * Ldft_land[0];      // Eq (A.3.1)
        Ldft[1] = omega * Ldft_sea[1] + (1 - omega) * Ldft_land[1];      // Eq (A.3.1)

        return Ldft;
    }

    double[] dl_se_ft_inner(double epsr, double sigma, double d, double hte, double hre, double adft, double f) {
        ////dl_se_ft_inner The inner routine of the first-term spherical diffraction loss
        //   This function computes the first-term part of Spherical-Earth diffraction
        //   loss exceeded for p// time for antenna heights
        //   as defined in Sec. A.3 of the ITU-R P.2001-3, equations (A3.2-A3.8)
        //
        //     Input parameters:
        //     epsr    -   Relative permittivity
        //     sigma   -   Conductivity (S/m)
        //     d       -   Great-circle path distance (km)
        //     hte     -   Effective height of interfering antenna (m)
        //     hre     -   Effective height of interfered-with antenna (m)
        //     adft    -   effective Earth radius (km)
        //     f       -   frequency (GHz)
        //
        //     Output parameters:
        //     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p// time
        //                implementing equations (30-37), Ldft(1) is for horizontal
        //                and Ldft(2) for the vertical polarization

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation (P.452)
        //     v1    06JUL16     Ivica Stevanovic, OFCOM         Modifications for P.1812
        //     v2    13JUL16     Ivica Stevanovic, OFCOM         Modifications for P.2001
        //     v3    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version
        ////

        // Normalized factor for surface admittance for horizontal (1) and vertical
        // (2) polarizations
        double[] K = new double[2];

        K[0] = 0.036 * Math.pow((adft * f), (-1.0 / 3.0)) * Math.pow((Math.pow((epsr - 1), 2) + Math.pow((18 * sigma / f), 2)), (-1.0 / 4.0));   // Eq (A.3.2a)

        K[1] = K[0] * Math.pow((Math.pow(epsr, 2) + Math.pow((18 * sigma / f), 2)), (1.0 / 2.0));       // Eq (A.3.2b)

        // Earth ground/polarization parameter
        double[] beta_dft = new double[2];
        double[] X = new double[2];
        double[] Yt = new double[2];
        double[] Yr = new double[2];

        for (int ii = 0; ii < 2; ii++) {
            beta_dft[ii] = (1 + 1.6 * Math.pow(K[ii], 2) + 0.67 * Math.pow(K[ii], 4)) / (1 + 4.5 * Math.pow(K[ii], 2) + 1.53 * Math.pow(K[ii], 4));  // Eq (A.3.3)

            // Normalized distance

            X[ii] = 21.88 * beta_dft[ii] * Math.pow((f / Math.pow(adft, 2)), (1.0 / 3.0)) * d;          // Eq (A.3.4)

            // Normalized transmitter and receiver heights

            Yt[ii] = 0.9575 * beta_dft[ii] * Math.pow((f * f / adft), (1.0 / 3.0)) * hte;       // Eq (A.3.5a)

            Yr[ii] = 0.9575 * beta_dft[ii] * Math.pow((f * f / adft), (1.0 / 3.0)) * hre;       // Eq (A.3.5b)
        }

        double[] Fx = new double[]{0, 0};
        double[] GYt = new double[]{0, 0};
        double[] GYr = new double[]{0, 0};
        double[] Bt = new double[]{0, 0};
        double[] Br = new double[]{0, 0};
        double[] Ldft = new double[]{0, 0};

        // Calculate the distance term given by:

        for (int ii = 0; ii < 2; ii++) {
            if (X[ii] >= 1.6) {
                Fx[ii] = 11 + 10 * Math.log10(X[ii]) - 17.6 * X[ii];
            } else {
                Fx[ii] = -20 * Math.log10(X[ii]) - 5.6488 * Math.pow((X[ii]), 1.425);     // Eq (A.3.6)
            }

            Bt[ii] = beta_dft[ii] * Yt[ii];             // Eq (A.3.7a)

            Br[ii] = beta_dft[ii] * Yr[ii];              // Eq (A.3.7a)

            if (Bt[ii] > 2) {
                GYt[ii] = 17.6 * Math.pow((Bt[ii] - 1.1), 0.5) - 5 * Math.log10(Bt[ii] - 1.1) - 8;
            } else {
                GYt[ii] = 20 * Math.log10(Bt[ii] + 0.1 * Math.pow(Bt[ii], 3));
            }

            if (Br[ii] > 2) {
                GYr[ii] = 17.6 * Math.pow((Br[ii] - 1.1), 0.5) - 5 * Math.log10(Br[ii] - 1.1) - 8;
            } else {
                GYr[ii] = 20 * Math.log10(Br[ii] + 0.1 * Math.pow(Br[ii], 3));
            }

            if (GYr[ii] < 2 + 20 * Math.log10(K[ii])) {
                GYr[ii] = 2 + 20 * Math.log10(K[ii]);
            }

            if (GYt[ii] < 2 + 20 * Math.log10(K[ii])) {
                GYt[ii] = 2 + 20 * Math.log10(K[ii]);
            }


            Ldft[ii] = -Fx[ii] - GYt[ii] - GYr[ii];         // Eq (A.3.8)

        }

        return Ldft;
    }

    public double[] dl_bull_actual(double[] d, double[] h, double hts, double hrs, double Cp, double f) {
        //dl_bull_actual Bullington part of the diffraction loss according to P.2001-3
        //   This function computes the Bullington part of the diffraction loss
        //   as defined in ITU-R P.2001-3 in Attachment A.4 (for the smooth profile)
        //
        //     Input parameters:
        //     d       -   Vector of distances di of the i-th profile point (km)
        //     h       -   Vector of heights hi of the i-th profile point (meters
        //                 above mean sea level)
        //                 Both vectors d and h contain n+1 profile points
        //     hts     -   Effective transmitter antenna height in meters above sea level (i=0)
        //     hrs     -   Effective receiver antenna height in meters above sea level (i=n)
        //     Cp      -   Effective Earth curvature
        //     f       -   Frequency (GHz)
        //
        //     Output parameters:
        //     Ldba   -   Bullington diffraction loss for a given actual path
        //     Ldbka  -   Knife-edge diffraction loss for Bullington point: actual path
        //     FlagLospa - 1 = LoS p// time for actual path, 0 = otherwise


        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation for P.452-16
        //     v1    06JUL16     Ivica Stevanovic, OFCOM         First implementation for P.1812-4
        //     v2    13JUL16     Ivica Stevanovic, OFCOM         First implementation for P.2001-3
        //     v3    13JUL17     Ivica Stevanovic, OFCOM         Initial Java implementation

        double Ldba = 0;
        double Ldbka = 0;
        int FlagLospa = 0;

        // Wavelength
        double c0 = 2.998e8;
        double lam = 1e-9 * c0 / f;

        // Complete path length
        int n = d.length;
        double dtot = d[n - 1] - d[0];

        // Find the intermediate profile point with the highest slope of the line
        // from the transmitter to the point

        int ii = 1;
        double Stim = ((h[ii] + 500 * Cp * d[ii] * (dtot - d[ii]) - hts) / d[ii]);           // Eq (A.4.1)

        for (ii = 2; ii < n - 1; ii++) {
            Stim = Math.max(Stim, ((h[ii] + 500 * Cp * d[ii] * (dtot - d[ii]) - hts) / d[ii]));           // Eq (A.4.1)
        }
        // Calculate the slope of the line from transmitter to receiver assuming a
        // LoS path

        double Str = (hrs - hts) / dtot;                                         // Eq (A.4.2)

        if (Stim < Str) { // Case 1, Path is LoS
            FlagLospa = 1;

            // Find the intermediate profile point with the highest diffraction
            // parameter nu:
            ii = 1;
            double numax = (h[ii] + 500 * Cp * d[ii] * (dtot - d[ii]) - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot) *
                    Math.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])));
            for (ii = 2; ii < n - 1; ii++) {

                double num = (h[ii] + 500 * Cp * d[ii] * (dtot - d[ii]) - (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot) *
                        Math.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])));   // Eq (A.4.3)
                if (num > numax) {
                    numax = num;
                }
            }
            Ldbka = dl_knife_edge(numax);                               // Eq (A.4.4)
        } else {
            FlagLospa = 0;
            // Path is NLOS

            // Find the intermediate profile point with the highest slope of the
            // line from the receiver to the point
            ii = 1;
            double Srim = (h[ii] + 500 * Cp * d[ii] * (dtot - d[ii]) - hrs) / (dtot - d[ii]);
            for (ii = 2; ii < n - 1; ii++) {
                double Sri = (h[ii] + 500 * Cp * d[ii] * (dtot - d[ii]) - hrs) / (dtot - d[ii]);   // Eq (A.4.5)
                if (Sri > Srim) {
                    Srim = Sri;
                }
            }
            // Calculate the distance of the Bullington point from the transmitter:

            double dbp = (hrs - hts + Srim * dtot) / (Stim + Srim);                // Eq (A.4.6)

            // Calculate the diffraction parameter, nub, for the Bullington point

            double nub = (hts + Stim * dbp - (hts * (dtot - dbp) + hrs * dbp) / dtot) *
                    Math.sqrt(0.002 * dtot / (lam * dbp * (dtot - dbp)));       // Eq (A.4.7)

            // The knife-edge loss for the Bullington point is given by

            Ldbka = dl_knife_edge(nub);                                 // Eq (A.4.8)

        }

        // For Luc calculated using either (A.4.4) or (A.4.8), Bullington diffraction loss
        // for the path is given by

        Ldba = Ldbka + (1 - Math.exp(-Ldbka / 6.0)) * (10 + 0.02 * dtot);            // Eq (A.4.9)

        double[] out = new double[3];

        out[0] = Ldba;
        out[1] = Ldbka;
        out[2] = (double) FlagLospa;

        return out;
    }

    double[] dl_bull_smooth(double[] d, double[] h, double htep, double hrep, double ap, double f) {
        //dl_bull_smooth Bullington part of the diffraction loss according to P.2001-3
        //   This function computes the Bullington part of the diffraction loss
        //   as defined in ITU-R P.2001-3 in Attachment A.5 (for a notional smooth profile)
        //
        //     Input parameters:
        //     d       -   Vector of distances di of the i-th profile point (km)
        //     h       -   Vector of heights hi of the i-th profile point (meters
        //                 above mean sea level)
        //                 Both vectors d and h contain n+1 profile points
        //     htep     -   Effective transmitter antenna height in meters above sea level (i=0)
        //     hrep     -   Effective receiver antenna height in meters above sea level (i=n)
        //     ap      -   Effective earth radius in kilometers
        //     f       -   frequency expressed in GHz
        //
        //     Output parameters:
        //     Ldbs   -   Bullington diffraction loss for a given smooth path
        //     Ldbks  -   Knife-edge diffraction loss for Bullington point: smooth path
        //     FlagLosps - 1 = LoS p// time for smooth path, 0 = otherwise

        //     Example:
        //     [Ldbs, Ldbks] = dl_bull_smooth(d, h, htep, hrep, ap, f)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation for P.452-16
        //     v1    06JUL16     Ivica Stevanovic, OFCOM         First implementation for P.1812-4
        //     v2    13JUL16     Ivica Stevanovic, OFCOM         First implementation for P.2001-3
        //     v3    29MAY17     Ivica Stevanovic, OFCOM         bug Jdbsk -> Ldbsk for NLOS path
        //                                                       as pointed out by Michael Rohner
        //     v4    13JUL17     Ivica Stevanovic, OFCOM         Initial Java implementation


        int FlagLosps = 1;
        double Ldbks = 0;

        // Wavelength
        double c0 = 2.998e8;
        double lam = 1e-9 * c0 / f;

        // Complete path length
        int n = d.length;
        double dtot = d[n - 1] - d[0];

        // Find the intermediate profile point with the highest slope of the line
        // from the transmitter to the point

        int ii = 1;
        double Stim = (500 * (dtot - d[ii]) / ap - htep / d[ii]);           // Eq (A.5.1)
        for (ii = 2; ii < n - 1; ii++) {
            double Sti = (500 * (dtot - d[ii]) / ap - htep / d[ii]);           // Eq (A.5.1)
            if (Sti > Stim) {
                Stim = Sti;
            }
        }
        // Calculate the slope of the line from transmitter to receiver assuming a
        // LoS path

        double Str = (hrep - htep) / dtot;                                         // Eq (A.5.2)

        if (Stim < Str) { // Case 1, Path is LoS
            FlagLosps = 1;
            // Find the intermediate profile point with the highest diffraction
            // parameter nu:
            ii = 1;
            double numax = (500 * d[ii] * (dtot - d[ii]) / ap - (htep * (dtot - d[ii]) + hrep * d[ii]) / dtot) *
                    Math.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])));

            for (ii = 2; ii < n - 1; ii++) {
                double num = (500 * d[ii] * (dtot - d[ii]) / ap - (htep * (dtot - d[ii]) + hrep * d[ii]) / dtot) *
                        Math.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])));     // Eq (A.5.3)
                if (num > numax) {
                    numax = num;
                }
            }
            Ldbks = dl_knife_edge(numax);                               // Eq (A.5.4)
        } else {
            FlagLosps = 0;
            // Path is NLOS

            // Find the intermediate profile point with the highest slope of the
            // line from the receiver to the point
            ii = 1;
            double Srim = (500 * d[ii] / ap - hrep / (dtot - d[ii]));     // Eq (A.5.5)

            for (ii = 2; ii < n - 1; ii++) {
                double Sri = (500 * d[ii] / ap - hrep / (dtot - d[ii]));
                if (Sri > Srim) {
                    Srim = Sri;
                }

            }

            // Calculate the distance of the Bullington point from the transmitter:

            double dbp = (hrep - htep + Srim * dtot) / (Stim + Srim);                // Eq (A.5.6)

            // Calculate the diffraction parameter, nub, for the Bullington point

            double nub = (htep + Stim * dbp - (htep * (dtot - dbp) + hrep * dbp) / dtot) *
                    Math.sqrt(0.002 * dtot / (lam * dbp * (dtot - dbp)));       // Eq (A.5.7)

            // The knife-edge loss for the Bullington point is given by

            Ldbks = dl_knife_edge(nub);                                 // Eq (A.5.8)

        }

        // For Ldbs calculated using either (A.5.4) or (A.5.8), Bullington diffraction loss
        // for the path is given by

        double Ldbs = Ldbks + (1 - Math.exp(-Ldbks / 6.0)) * (10 + 0.02 * dtot);            // Eq (A.5.9)

        double[] out = new double[3];
        out[0] = Ldbs;
        out[1] = Ldbks;
        out[2] = (double) FlagLosps;

        return out;
    }


    private double dl_knife_edge(double nu) {
        //dl_knife_edge Knife edge diffraction loss
        //   This function computes knife-edge diffraction loss in dB
        //   as defined in ITU-R P.2001-3 Section 3.12
        //
        //     Input parameters:
        //     nu      -   dimensionless parameter
        //
        //     Output parameters:
        //     J       -   Knife-edge diffraction loss (dB)
        //
        //     Example J = dl_knife_edge(nu)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    14JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        double J = 0;

        if (nu > -0.78) {
            J = 6.9 + 20 * Math.log10(Math.sqrt(Math.pow((nu - 0.1), 2) + 1) + nu - 0.1);
        }

        return J;
    }


    public double multi_path_activity(double f, double dt, double hts, double hrs, double dlt, double dlr, double hlt, double hlr, double hlo, double thetat, double thetar, double epsp, double Nd65m1, double phimn, int FlagLos50) {
        //multi_path_activity Multipath fading calculation
        //   This function computes the the notional zero-fade annual percentage time
        //   for the whole path as defined in ITU-R P.2001-3 in Attachment B.2
        //
        //     Input parameters:
        //     f        -   Frequency (GHz)
        //     hts      -   Transmitter antenna height in meters above sea level (i=0)
        //     hrs      -   Receiver antenna height in meters above sea level (i=n)
        //     dlt      -   Tx antenna horizon distance (km)
        //     dlr      -   Rx antenna horizon distance (km)
        //     hlt      -   Profile height at transmitter horizon (m)
        //     hlr      -   Profile height at receiver horizon (m)
        //     hlo      -   Lower antenna height (m)
        //     thetat   -   Horizon elevation angles relative to the local horizontal as viewed from Tx
        //     thetar   -   Horizon elevation angles relative to the local horizontal as viewed from Rx
        //     phimn    -   mid-point latitude
        //     dca      -   path distance (km)
        //     epsp     -   Positive value of path inclination (mrad)
        //     Nd65m1   -   Refractivity gradient in the lowest 65 m of the atmosphere exceeded for 1// of an average year
        //     phimn    -   Path midpoint latitude
        //     FlagLos50-   1 = Line-of-sight 50// time, 0 = otherwise
        //
        //     Output parameters:
        //     Q0ca   -   Notional zero-fade annual percantage time for the whole path
        //
        //
        //
        //     Example:
        //     Q0ca = multi_path_activity(f, dt, hts, hrs, dlt, dlr, hlt, hlr, hlo, thetat, thetar, epsp, Nd65m1, phimn)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version
        //// B.2 Characterize multi-path activity

        double Q0ca;

        // Factor representing the statistics of radio-refractivity lapse time

        double K = Math.pow(10, (-(4.6 + 0.0027 * Nd65m1)));

        if (FlagLos50 == 1) {
            // Calculate the notional zero fade annual percentage using (B.2.2)

            double dca = dt;
            double epsca = epsp;
            double hca = hlo;

            Q0ca = zero_fade_annual_time(dca, epsca, hca, f, K, phimn);


        } else {

            // Calculate the notoinal zero-fade annual percentage time at the
            // transmitter end Q0cat using (B.2.3)

            double dca = dlt;
            double epsca = Math.abs(thetat);
            double hca = Math.min(hts, hlt);

            double Q0cat = zero_fade_annual_time(dca, epsca, hca, f, K, phimn);

            // Calculate th enotional zero-fade annual percentage time at the
            // receiver end Q0car using (B.2.4)

            dca = dlr;
            epsca = Math.abs(thetar);
            hca = Math.min(hrs, hlr);

            double Q0car = zero_fade_annual_time(dca, epsca, hca, f, K, phimn);

            Q0ca = Math.max(Q0cat, Q0car);   // Eq (B.2.5)
        }

        return Q0ca;
    }


    public double zero_fade_annual_time(double dca, double epsca, double hca, double f, double K, double phimn) {
        //zero_fade_annual_time Calculate the notional zero-fade annual percentage time
        //   This function computes the the notional zero-fade annual percentage time
        //   as defined in ITU-R P.2001-3 in Attachment B.3
        //
        //     Input parameters:
        //     dca      -   path distance (km)
        //     epsca    -   Positive value of path inclination (mrad)
        //     hca      -   Antenna height in meters above sea level
        //     f        -   Frequency (GHz)
        //     K        -   Factor representing the statistics of radio-refractivity
        //                  lapse rate for the midpoint of the path
        //     phimn    -   mid-point latitude
        //     Output parameters:
        //     Q0ca   -   Notional zero-fade annual percantage time
        //
        //
        //     Example:
        //     Q0ca = zero_fade_annual_time(dca, epsca, hca, f, K, phimn)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version/
        //     v1    13JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        double Cg = 0;

        // Notional zero-fade worst-month percentage time (B.3.1)

        double qw = K * Math.pow(dca, 3.1) * Math.pow((1 + epsca), (-1.29)) * Math.pow(f, 0.8) * Math.pow(10, (-0.00089 * hca));

        // Calculate the logarithmic climatic conversion factor (B.3.2)

        if (Math.abs(phimn) <= 45) {

            Cg = 10.5 - 5.6 * Math.log10(1.1 + Math.pow((Math.abs(cosd(2 * phimn))), 0.7)) - 2.7 * Math.log10(dca) + 1.7 * Math.log10(1 + epsca);

        } else {

            Cg = 10.5 - 5.6 * Math.log10(1.1 - Math.pow((Math.abs(cosd(2 * phimn))), 0.7)) - 2.7 * Math.log10(dca) + 1.7 * Math.log10(1 + epsca);

        }

        if (Cg > 10.8) {
            Cg = 10.8;
        }

        // Notional zero-fade annual percentage time (B.3.3)

        double Q0ca = qw * Math.pow(10, (-0.1 * Cg));

        return Q0ca;
    }


    public double[] precipitation_fade_initial(P2001DigitalMaps maps, double f, double q, double phi_n, double phi_e, double h_rainlo, double h_rainhi, double d_rain, int pol_hv) {
        //precipitation_fade_initial Preliminary calculation of precipitation fading
        //   This function computes the preliminary parameters necessary for precipitation
        //   fading as defined in ITU-R P.2001-3 in Attachment C.2
        //
        //     Input parameters:
        //     maps     -   Object containing digital maps
        //     f        -   Frequency (GHz)
        //     q        -   Percentage of average year for which predicted basic
        //                  loss is exceeded (100-p)
        //     phi_n    -   Latitude for obtaining rain climatic parameters (deg)
        //     phi_e    -   Longitude for obtaining rain climatic parameters (deg)
        //     h_rainlo -   Lower height of the end of the path for a precipitation calculation (m)
        //     h_rainlo -   Higher height of the end of the path for a precipitation calculation (m)
        //     d_rain   -   Length of the path for rain calculation (km)
        //     pol_hv   -   0 = horizontal, 1 = vertical polarization
        //
        //
        //     Output parameters:
        //     a, b, c  -   Parameters defining cumulative distribution of rain rate
        //     dr       -   Limited path length for precipitation calculations
        //     Q0ra     -   Percentage of an average year in which rain occurs
        //     Fwvr     -   Factor used to estimate the effect of additional water vapour under rainy conditions
        //     kmod     -   Modified regression coefficients
        //     alpha_mod-   Modified regression coefficients
        //     ngp      -   size of vector G (i.e., P)
        //     G        -   Vector of attenuation multipliers
        //     P        -   Vector of probabilities
        //     flagrain -   1 = "rain" path, 0 = "non-rain" path
        //
        //     Example:
        //     [a, b, c, dr, Q0ra, Fwvr, kmod, alpha_mod, ngp, G, P, flagrain] = precipitation_fade_initial(f, q, phi_n, phi_e, h_rainlo, h_rainhi, d_rain, pol_hv)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    15JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUN17     Ivica Stevanovic, OFCOM         replaced load calls to increase computational speed
        //     v2    19JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        //// C.2 Precipitation fading: Preliminary calculations

        // Obtain Pr6 for phi_n, phi_e from the data file "Esarain_Pr6_v5.txt"
        // as a bilinear interpolation

        //Esarain_Pr6 = load('DigitalMaps/Esarain_Pr6_v5.txt');
        double Pr6 = maps.GetPr6(phi_e, phi_n);


        // Obtain Mt for phi_n, phi_e from the data file "Esarain_Mt_v5.txt"
        // as a bilinear interpolation

        //Esarain_Mt = load('DigitalMaps/Esarain_Mt_v5.txt');
        double Mt = maps.GetMt(phi_e, phi_n);

        // Obtain beta_rain for phi_n, phi_e from the data file "Esarain_Beta_v5.txt"
        // as a bilinear interpolation

        //Esarain_Beta = load('DigitalMaps/Esarain_Beta_v5.txt');
        double beta_rain = maps.GetBetaRain(phi_e, phi_n);

        // Obtain h0 for phi_n, phi_e from the data file "h0.txt"
        // as a bilinear interpolatio

        double h0 = maps.GetH0(phi_e, phi_n);

        // Calculate mean rain height hr (C.2.1)

        double hR = 360 + 1000 * h0;

        // calculate the highest rain height hRtop (C.2.2)

        double hRtop = hR + 2400;

        double flagrain; // no rain path
        double Q0ra;
        double Fwvr;
        double a;
        double b;
        double c;
        double dr;
        double kmod;
        double k;
        double alpha_mod;
        double alpha;
        double[] G;
        double[] P;
        int ngp;


        if (Pr6 == 0 || h_rainlo >= hRtop) {
            flagrain = 0; // no rain path
            Q0ra = 0;
            Fwvr = 0;
            a = 0;
            b = 0;
            c = 0;
            dr = 0;
            kmod = 0;
            alpha_mod = 0;
            G = new double[1];
            P = new double[1];
            G[0] = 0;
            P[0] = 0;
            ngp = 1;
        } else { // the path is classified as rain

            flagrain = 1;

            // Values from table C.2.1

            double[] H = {-2400.0,
                    -2300.0,
                    -2200.0,
                    -2100.0,
                    -2000.0,
                    -1900.0,
                    -1800.0,
                    -1700.0,
                    -1600.0,
                    -1500.0,
                    -1400.0,
                    -1300.0,
                    -1200.0,
                    -1100.0,
                    -1000.0,
                    -900.0,
                    -800.0,
                    -700.0,
                    -600.0,
                    -500.0,
                    -400.0,
                    -300.0,
                    -200.0,
                    -100.0,
                    0.0,
                    100.0,
                    200.0,
                    300.0,
                    400.0,
                    500.0,
                    600.0,
                    700.0,
                    800.0,
                    900.0,
                    1000.0,
                    1100.0,
                    1200.0,
                    1300.0,
                    1400.0,
                    1500.0,
                    1600.0,
                    1700.0,
                    1800.0,
                    1900.0,
                    2000.0,
                    2100.0,
                    2200.0,
                    2300.0,
                    2400.0
            };

            double[] Pi = {0.000555,
                    0.000802,
                    0.001139,
                    0.001594,
                    0.002196,
                    0.002978,
                    0.003976,
                    0.005227,
                    0.006764,
                    0.008617,
                    0.010808,
                    0.013346,
                    0.016225,
                    0.019419,
                    0.022881,
                    0.026542,
                    0.030312,
                    0.034081,
                    0.037724,
                    0.041110,
                    0.044104,
                    0.046583,
                    0.048439,
                    0.049589,
                    0.049978,
                    0.049589,
                    0.048439,
                    0.046583,
                    0.044104,
                    0.041110,
                    0.037724,
                    0.034081,
                    0.030312,
                    0.026542,
                    0.022881,
                    0.019419,
                    0.016225,
                    0.013346,
                    0.010808,
                    0.008617,
                    0.006764,
                    0.005227,
                    0.003976,
                    0.002978,
                    0.002196,
                    0.001594,
                    0.001139,
                    0.000802,
                    0.000555
            };

            // Calculate two intermediate parameters (C.2.3)

            double Mc = beta_rain * Mt;
            double Ms = (1 - beta_rain) * Mt;

            // Calculate the percentage of an average year in which rain occurs (C.2.4)

            Q0ra = Pr6 * (1 - Math.exp(-0.0079 * Ms / Pr6));

            // Calculate the parameters defining the cumulative distribution of rain
            // rate (C.2.5)

            double a1 = 1.09;
            double b1 = (Mc + Ms) / (21797 * Q0ra);
            double c1 = 26.02 * b1;

            a = a1;
            b = b1;
            c = c1;

            // Calculate the percentage time approximating to the transition between
            // the straight and curved sections of the rain-rate cumulative
            // distribution when plotted ... (C.2.6)

            double Qtran = Q0ra * Math.exp(a1 * (2 * b1 - c1) / (c1 * c1));


            // Path inclination angle (C.2.7)

            double eps_rain = 0.001 * (h_rainhi - h_rainlo) / d_rain; // radians

            // Use the method given in Recommendation ITU-R P.838 to calculate rain
            // regression coefficients k and alpha for the frequency, polarization
            // and path inclination

            if (f < 1) { // compute theregression coefficient for 1 GHz

                double[] pp = p838(1, eps_rain, pol_hv);
                double k1GHz = pp[0];
                double alpha1GHz = pp[1];
                k = f * k1GHz;        // Eq (C.2.8a)
                alpha = alpha1GHz;  // Eq (C.2.8b)
            } else {

                double[] pp = p838(f, eps_rain, pol_hv);
                k = pp[0];
                alpha = pp[1];
            }

            // Limit the path length for precipitation (C.2.9)

            dr = Math.min(d_rain, 300.0);
            double drmin = Math.max(dr, 1.0);

            // Calculate modified regression coefficients (C.2.10)

            kmod = Math.pow(1.763, alpha) * k * (0.6546 * Math.exp(-0.009516 * drmin) + 0.3499 * Math.exp(-0.001182 * drmin));
            alpha_mod = (0.753 + 0.197 / drmin) * alpha + 0.1572 * Math.exp(-0.02268 * drmin) - 0.1594 * Math.exp(-0.0003617 * drmin);

            // Initialize and allocate the arrays for attenuation multiplier and
            // probability of a particular case (with a maximum dimension 49 as in
            // table C.2.1

            int nn = H.length;
            double[] Gm = new double[nn];
            double[] Pm = new double[nn];

            // Initialize Gm(1)=1, set m=1

            Gm[0] = 1;

            int m = 0;

            for (int n = 0; n < 49; n++) {
                // For each line of Table C.2.1 for n from 1 to 49 do the following

                // a) Calculate rain height given by  (C.2.11)
                double hT = hR + H[n];

                // b) If h_rainlo >= hT repeat from a) for the next avalue of n,
                // otherwise continue from c

                if (h_rainlo >= hT) {

                    continue; // repeat from a) for the next value of n
                }


                if (h_rainhi > hT - 1200) {
                    // c.i) use the method in Attachment C.5 to set Gm to the
                    // path-averaged multiplier for this path geometry relative to
                    // the melting layer

                    Gm[m] = path_averaged_multiplier(h_rainlo, h_rainhi, hT);

                    // c.ii) set Pm = Pi(n) from Table C.2.1
                    Pm[m] = Pi[n];

                    // c.iii) if n < 49 add 1 to array index m

                    if (n < 49 - 1) {
                        m = m + 1;
                    }

                    // c.iv) repeat fom a) for the next value of n

                    //continue;

                } else {

                    // d) Accumulate Pi(n) from table C.2.1 into Pm, set Gm = 1 and
                    // repeat from a) for the next value of n

                    Pm[m] = Pm[m] + Pi[n];
                    Gm[m] = 1;
                }

            } // for loop

            // Set the number of values in arrays Gm and Pm according to (C.2.12)

            int Mlen = m + 1;

            G = new double[Mlen];
            P = new double[Mlen];
            for (int n = 0; n < Mlen; n++) {
                G[n] = Gm[n];
                P[n] = Pm[n];
            }
            ngp = Mlen;

            // Calculate a factor used to estimate the effect of additional water
            // vapour under rainy conditions (C.2.13), (C.2.14)

            double Rwvr = 6 * (Math.log10(Q0ra / q) / Math.log10(Q0ra / Qtran)) - 3;

            double GP = 0;
            for (int n = 0; n < G.length; n++) {
                GP = GP + G[n] * P[n];
            }
            Fwvr = 0.5 * (1 + Math.tanh(Rwvr)) * GP;

        }


        // [a, b, c, dr, Q0ra, Fwvr, kmod, alpha_mod, ngp, G, P, flagrain]
        double[] out = new double[2 * ngp + 10];
        out[0] = a;
        out[1] = b;
        out[2] = c;
        out[3] = dr;
        out[4] = Q0ra;
        out[5] = Fwvr;
        out[6] = kmod;
        out[7] = alpha_mod;
        out[8] = ngp;
        for (int i = 0; i < ngp; i++) {
            out[9 + i] = G[i];
            out[9 + ngp + i] = P[i];
        }
        out[9 + 2 * ngp] = flagrain;
        return out;
    }


    public double[] p838(double f, double theta, int pol) {
        //p838 Recommendation ITU-R P.838-3
        //   This function computes the rain regression coefficients k and alpha for
        //   a given frequency, path inclination and polarization according to ITU-R
        //   Recommendation P.838-3
        //
        //     Input parameters:
        //     f        -   Frequency (GHz)
        //     theta    -   Path inclination (radians)
        //     pol      -   Polarization 0 = horizontal, 1 = vertical
        //
        //     Output parameters:
        //     k, alpha -   Rain regression coefficients
        //
        //     Example:
        //     [ k, alpha ] = p838( f, theta, pol )

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    17JUL17     Ivica Stevanovic, OFCOM         Initial Java implementation

        double tau = 0;

        if (pol == 0) { // horizontal polarization

            tau = 0;

        } else { // vertical polarization

            tau = Math.PI / 2.0;

        }

        // Coefficients for kH
        double[] aj_kh = {-5.33980,
                0.35351,
                0.23789,
                0.94158
        };

        double[] bj_kh = {-0.10008,
                -1.26970,
                -0.86036,
                -0.64552
        };

        double[] cj_kh = {1.13098,
                0.45400,
                0.15354,
                0.16817
        };


        double m_kh = -0.18961;
        double c_kh = 0.71147;

        // Coefficients for kV

        double[] aj_kv = {
                -3.80595,
                -3.44965,
                -0.39902,
                0.50167
        };

        double[] bj_kv = {
                0.56934,
                -0.22911,
                0.73042,
                1.07319
        };
        double[] cj_kv = {
                0.81061,
                0.51059,
                0.11899,
                0.27195,
        };
        double m_kv = -0.16398;
        double c_kv = 0.63297;

        // Coefficients for aH


        double[] aj_ah = {
                -0.14318,
                0.29591,
                0.32177,
                -5.37610,
                16.1721
        };

        double[] bj_ah = {
                1.82442,
                0.77564,
                0.63773,
                -0.96230,
                -3.29980
        };

        double[] cj_ah = {
                -0.55187,
                0.19822,
                0.13164,
                1.47828,
                3.43990
        };
        double m_ah = 0.67849;
        double c_ah = -1.95537;

        // Coefficients for aV

        double[] aj_av = {
                -0.07771,
                0.56727,
                -0.20238,
                -48.2991,
                48.5833
        };

        double[] bj_av = {
                2.33840,
                0.95545,
                1.14520,
                0.791669,
                0.791459
        };

        double[] cj_av = {
                -0.76284,
                0.54039,
                0.26809,
                0.116226,
                0.116479
        };

        double m_av = -0.053739;
        double c_av = 0.83433;
        double logkh = 0;

        for (int j = 0; j < aj_kh.length; j++) {
            logkh = logkh + aj_kh[j] * Math.exp(-Math.pow(((Math.log10(f) - bj_kh[j]) / cj_kh[j]), 2.0));
        }
        logkh = logkh + m_kh * Math.log10(f) + c_kh;

        double kh = Math.pow(10.0, (logkh));

        double logkv = 0.0;

        for (int j = 0; j < aj_kv.length; j++) {
            logkv = logkv + aj_kv[j] * Math.exp(-Math.pow((Math.log10(f) - bj_kv[j]) / cj_kv[j], 2.0));
        }
        logkv = logkv + m_kv * Math.log10(f) + c_kv;

        double kv = Math.pow(10.0, (logkv));

        double ah = 0.0;

        double av = 0.0;

        for (int j = 0; j < aj_ah.length; j++) {
            ah = ah + aj_ah[j] * Math.exp(-Math.pow(((Math.log10(f) - bj_ah[j]) / cj_ah[j]), 2.0));

            av = av + aj_av[j] * Math.exp(-Math.pow(((Math.log10(f) - bj_av[j]) / cj_av[j]), 2.0));
        }

        ah = ah + m_ah * Math.log10(f) + c_ah;

        av = av + m_av * Math.log10(f) + c_av;

        double k = (kh + kv + (kh - kv) * Math.pow((Math.cos(theta)), 2.0) * Math.cos(2.0 * tau)) / 2.0;

        double alpha = (kh * ah + kv * av + (kh * ah - kv * av) * Math.pow((Math.cos(theta)), 2.0) * Math.cos(2.0 * tau)) / (2.0 * k);

        double[] out = new double[2];
        // [ k, alpha ]
        out[0] = k;
        out[1] = alpha;

        return out;
    }


    public double path_averaged_multiplier(double hlo, double hhi, double hT) {
        //path_averaged_multiplier Models path-averaged multiplier
        //   This function computes the path-averaged multiplier
        //   according to ITU-R Recommendation P.2001-3 Attachment C.5
        //
        //     Input parameters:
        //     hlo, hhi  -   heights of the loewr and higher antennas (m)
        //     hT        -   rain height
        //
        //     Output parameters:
        //     G         -   weighted average of the multiplier Gamma (multi_layer.m)
        //
        //     Example:
        //     G = path_averaged_multiplier( hlo, hhi, hT )

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    17JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        //// C.5 Path-averaged multiplier

        // Calculate the slices in which the two antennas lie (C.5.1)

        double G;

        int slo = (int) (1 + Math.floor((hT - hlo) / 100));

        int shi = (int) (1 + Math.floor((hT - hhi) / 100));

        if (slo < 1) { // path wholly above the melting layer
            G = 0;
            return G;
        }

        if (shi > 12) { // path is wholy at or below the lower edge of the melting layer
            G = 1;
            return G;
        }

        if (slo == shi) { // both antennas in the same melting-layer slice (C.5.2)

            G = multi_layer(0.5 * (hlo + hhi) - hT);

            return G;
        }

        // Initialize G for use as an accumulator (C.5.3)

        G = 0;

        // Calculate the required range of slice indices (C.5.4)

        int sfirst = Math.max(shi, 1);

        int slast = Math.min(slo, 12);

        double delh = 0;
        double Q = 0;

        for (int s = sfirst; s <= slast; s++) {

            if (shi < s && s < slo) {
                // In this case the slice is fully traversed by a section of the
                // path (C.5.5)

                delh = 100 * (0.5 - s);

                Q = 100 / (hhi - hlo);

            } else if (s == slo) {
                // In this case the slice contains the lower antenna at hlo (C.5.6)

                delh = 0.5 * (hlo - hT - 100 * (s - 1));
                Q = (hT - 100 * (s - 1) - hlo) / (hhi - hlo);

            } else if (s == shi) {
                // In this case the slice contains the higher antenna at hhi (C.5.7)

                delh = 0.5 * (hhi - hT - 100 * s);

                Q = (hhi - (hT - 100 * s)) / (hhi - hlo);

            }

            // For delh calculated under one of the preceeding three conditions,
            // calculate the corresponding multiplier (C.5.8)

            double Gamma_slice = multi_layer(delh);

            // Accumulate the multiplier (C.5.9)

            G = G + Q * Gamma_slice;
        }

        if (slo > 12) { // lower antenna is below the melting layer

            //The fraction of the path below the layer (C.5.10)
            Q = (hT - 1200 - hlo) / (hhi - hlo);

            // Since the multiplier is 1 below the layer, G should be increased
            // according to (C.5.11)

            G = G + Q;

        }

        return G;
    }


    double multi_layer(double delh) {
        //multi_layer Models the changes in specific attenuation within the melting layer
        //   This function computes the changes in specific attenuation at different
        //   heights within the melting layer according to ITU-R Recommendation
        //   P.2001-3 Attachment C.4
        //
        //     Input parameters:
        //     delh     -   difference between a given height h and rain height hT (m)
        //
        //     Output parameters:
        //     Gamma    -   attenuation multiplier
        //
        //     Example:
        //     Gamma = multi_layer( delh )

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    17JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // C.4 Melting-layer model

        double Gamma;

        if (delh > 0) {

            Gamma = 0;

        } else if (delh < -1200) {

            Gamma = 1;

        } else {

            Gamma = 4 * Math.pow((1 - Math.exp(delh / 70)), 2);

            Gamma = Gamma / (1 + Math.pow((1 - Math.exp(-Math.pow((delh / 600), 2))), 2) * (Gamma - 1));

        }

        return Gamma;
    }


    public double Aiter(double q, double Q0ca, double Q0ra, int flagtropo, double a, double b, double c, double dr, double kmod, double alpha_mod, double[] Gm, double[] Pm, int flagrain) {
        //Aiter Inverse cumulative distribution function of a propagation model
        //   This function computes the inverse cumulative distribution function of a
        //   propagation model as defined in ITU-R P.2001-3 in Attachment I
        //
        //     Input parameters:
        //     q        -   Percentage of average year for which predicted
        //                  transmission loss is exceeded
        //     Q0ca     -   Notional zero-fade annual percentage time
        //     Q0ra     -   Percentage of an average year in which rain occurs
        //     flagtropo-   0 = surface path, 1 = tropospheric path
        //     a, b, c  -   Parameters defining cumulative distribution of rain rate
        //     dr       -   Limited path length for precipitation calculations
        //     kmod     -   Modified regression coefficients
        //     alpha_mod-   Modified regression coefficients
        //     Gm       -   Vector of attenuation multipliers
        //     Pm       -   Vector of probabilities
        //     flagrain -   1 = "rain" path, 0 = "non-rain" path
        //
        //     Output parameters:
        //     AiterQ   -   Attenuation level of a propoagation mechanisms exceeded
        //                  for q% time
        //
        //
        //
        //     Example:
        //     AiterQ = Aiter(q, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    15JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    17JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // Attachment I. Iterative procedure to invert cumulative distribution function

        // Set initial values of the high and low searc hlimits for attenuation and
        // the attenuation step size

        double Ainit = 10;

        double Ahigh = Ainit / 2;                                                          // Eq (I.2.1)

        double Alow = -Ainit / 2;                                                          // Eq (I.2.2)

        double Astep = Ainit;                                                            // Eq (I.2.3)

        // Initialize the percentage times attenuations Ahigh and Alow are exceeded
        // (I.2.4)

        double qhigh = Qiter(Ahigh, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);

        double qlow = Qiter(Alow, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);

        // Stage 1

        int count = 0;

        if (q < qhigh || q > qlow) {

            while (count < 11) {
                count = count + 1;

                if (q < qhigh) {

                    Alow = Ahigh;
                    qlow = qhigh;
                    Astep = 2 * Astep;
                    Ahigh = Ahigh + Astep;
                    qhigh = Qiter(Ahigh, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);
                    continue; // Loop back to the start of search range iteration and repeat from there
                }

                if (q > qlow) {

                    Ahigh = Alow;
                    qhigh = qlow;
                    Astep = 2 * Astep;
                    Alow = Alow - Astep;
                    qlow = Qiter(Alow, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);
                    continue; // Loop back to the start of search range iteration and repeat from there
                }

            } // Initial search range iteration

        }  // only if q < qhigh and q > qlow


        // Stage 2: Binary search

        // Evaluate Atry (I.2.5)

        double Atry = 0.5 * (Alow + Ahigh);

        // Start of binary search iteration
        // Set the binary search accuracy

        double Aacc = 0.01;

        int Niter = (int) Math.ceil(3.32 * Math.log10(Astep / Aacc));

        count = 0;
        double qtry;

        while (count <= Niter) {
            count = count + 1;

            // Calculate the percentage time attenuation Atry is exceeded (I.2.6)

            qtry = Qiter(Atry, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);

            if (qtry < q) {
                Ahigh = Atry;
            } else {
                Alow = Atry;
            }

            Atry = 0.5 * (Alow + Ahigh);

        }  // Loop back to the start of binary search iteration and repeat from there


        return Atry;
    }

    public double Qiter(double Afade, double Q0ca, double Q0ra, int flagtropo, double a, double b, double c, double dr, double kmod, double alpha_mod, double[] Gm, double[] Pm, int flagrain) {
        //Qiter Cumulative distribution function of a propagation model
        //   This function computes the cumulative distribution function of a
        //   propagation model as defined in ITU-R P.2001-3 in Sections 4.1 and 4.3
        //
        //     Input parameters:
        //     Afade    -   Clear air fade (A>0) or enhancement (A<0)
        //     Q0ca     -   Notional zero-fade annual percantage time
        //     Q0ra     -   Percentage of an average year in which rain occurs
        //     flagtropo-   0 = surface path, 1 = tropospheric path
        //     a, b, c  -   Parameters defining cumulative distribution of rain rate
        //     dr       -   Limited path length for precipitation calculations
        //     kmod     -   Modified regression coefficients
        //     alpha_mod-   Modified regression coefficients
        //     Gm       -   Vector of attenuation multipliers
        //     Pm       -   Vector of probabilities
        //     flagrain -   1 = "rain" path, 0 = "non-rain" path
        //
        //     Output parameters:
        //     QiterA   -   Cumulative distribution function of fade A
        //
        //
        //
        //     Example:
        //     QiterA = Qiter(Afade, Q0ca, Q0ra, flagtropo, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    17JUL17     Ivica Stevanovic, OFCOM         Initial Java version


        // Compute Qrain(Afade) as defined in Attachment C.3
        double QrainA = precipitation_fade(Afade, a, b, c, dr, kmod, alpha_mod, Gm, Pm, flagrain);

        // Compute Qcaf(Afade) as defined in Attachments B.4/5
        double QcafA;
        if (flagtropo == 0) {
            QcafA = clear_air_fade_surface(Afade, Q0ca);
        } else {
            QcafA = clear_air_fade_tropo(Afade);
        }

        // Function QiterA is defined for combined clear-air/precipitation fading
        // (4.1.3), (4.3.5)

        return QrainA * (Q0ra / 100) + QcafA * (1 - Q0ra / 100);


    }

    public double precipitation_fade(double Afade, double a, double b, double c, double dr, double kmod, double alpha_mod, double[] Gm, double[] Pm, int flagrain) {
        //precipitation_fade Percentage time for which attenuation is exceeded
        //   This function computes the percentage time during which it is raining
        //   for which a given attenuation Afade is exceeded as defined in ITU-R
        //   P.2001-3 in Attachment C.3
        //
        //     Input parameters:
        //     Afade    -   Attenuation (dB)
        //     a, b, c  -   Parameters defining cumulative distribution of rain rate
        //     dr       -   Limited path length for precipitation calculations
        //     kmod     -   Modified regression coefficients
        //     alpha_mod-   Modified regression coefficients
        //     Gm       -   Vector of attenuation multipliers
        //     Pm       -   Vector of probabilities
        //     flagrain -   1 = "rain" path, 0 = "non-rain" path
        //
        //     Output parameters
        //     QrainA   -  percentage time during which it is raining for which a
        //                 given attenuation Afade is exceeded
        //     Example:
        //     QrainA = precipitation_fade(Afade, a, b, c, dr, Q0ra, Fwvr, kmod, alpha_mod, Gm, Pm, flagrain)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    17JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // C.3 Percentage time a given precipitation fade level is exceeded

        double QrainA = 0;

        if (Afade < 0) {

            QrainA = 100;                                                        // Eq (C.3.1a)

        } else {

            if (flagrain == 0) {  // non-rain

                QrainA = 0;                                                       // Eq (C.3.1b)

            } else {    // rain

                double drlim = Math.max(dr, 0.001); // Eq (C.3.1e)
                int n = Gm.length;
                double Rm;
                QrainA = 0;
                for (int i = 0; i < n; i++) {

                    Rm = Math.pow((Afade / (Gm[i] * drlim * kmod)), (1 / alpha_mod));            // Eq (C.3.1d)

                    QrainA = QrainA + 100 * (Pm[i] * Math.exp(-a * Rm * (b * Rm + 1) / (c * Rm + 1)));        // Eq (C.3.1c)

                }
            }
        }
        return QrainA;
    }

    public double clear_air_fade_surface(double A, double Q0ca) {
        //Percentage time a given clear-air fade level is exceeded on a surface path
        //   This function computes the percentage of the non-rain time a given fade
        //   in dB below the median signal level is exceeded on a surface path
        //   as defined in ITU-R P.2001-3 in Attachment B.4
        //
        //     Input parameters:
        //     A        -   Clear air fade (A>0) or enhancement (A<0)
        //     Q0ca     -   Notional zero-fade annual percantage time
        //     Output parameters:
        //     Qcaf     -   Percentage time a given clear-air fade level is exceeded
        //                  on a surface path
        //
        //
        //
        //     Example:
        //     Qcaf = clear_air_fade_surface(A, Q0ca)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // B.4 Percentage time a given clear-air fade level is exceeded on a surface path

        double Qcaf;

        if (A >= 0) {


            double qt = 3.576 - 1.955 * Math.log10(Q0ca);     // Eq (B.4.1b)

            double qa = 2 + (1 + 0.3 * Math.pow(10, (-0.05 * A))) * Math.pow(10, (-0.016 * A)) * (qt + 4.3 * (Math.pow(10, (-0.05 * A)) + A / 800));    // Eq (B.4.1a)

            Qcaf = 100 * (1 - Math.exp(-Math.pow(10, (-0.05 * qa * A)) * Math.log(2)));   // Eq (B.4.1)

        } else {

            double qs = -4.05 - 2.35 * Math.log10(Q0ca);     // Eq (B.4.2b)

            double qe = 8 + (1 + 0.3 * Math.pow(10, (0.05 * A))) * Math.pow(10, (0.035 * A)) * (qs + 12 * (Math.pow(10, (0.05 * A)) - A / 800));    // Eq (B.4.2a)

            Qcaf = 100 * (Math.exp(-Math.pow(10, (0.05 * qe * A)) * Math.log(2)));   // Eq (B.4.2)
        }

        return Qcaf;
    }

    public double clear_air_fade_tropo(double A) {
        //Percentage time a given clear-air fade level is exceeded on a troposcatter path
        //   This function computes the percentage of the non-rain time a given fade
        //   in dB below the median signal level is exceeded on a troposcatter path
        //   as defined in ITU-R P.2001-3 in Attachment B.5
        //
        //     Input parameters:
        //     A        -   Clear air fade (A>0) or enhancement (A<0)
        //
        //     Output parameters:
        //     Qcaftropo-   Percentage time a given clear-air fade level is exceeded
        //                  on a troposcatter path
        //
        //
        //
        //     Example:
        //     Qcaftropo = clear_air_fade_tropo(A)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    14JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // B.5 Percentage time a given clear-air fade level is exceeded on a troposcatter path
        double Qcaftropo;
        if (A < 0) {

            Qcaftropo = 100;   // Eq (B.5.1a)

        } else {

            Qcaftropo = 0;     // Eq (B.5.1b)

        }

        return Qcaftropo;
    }

    private double tl_free_space(double f, double d) {
        //tl_free_space Free-space basic transmission loss
        //   This function computes free-space basic transmission loss in dB
        //   as defined in ITU-R P.2001-3 Section 3.11
        //
        //     Input parameters:
        //     f       -   Frequency (GHz)
        //     d       -   Distance (km)
        //
        //     Output parameters:
        //     Lbfs    -   Free-space basic transmission loss (dB)

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        return 92.44 + 20 * Math.log10(f) + 20 * Math.log10(d);

    }

    public double[] tl_anomalous_reflection(double f, double[] d, int[] z, double hts, double hrs, double htea, double hrea, double hm, double thetat, double thetar, double dlt, double dlr, double phimn, double omega, double ae, double p, double q) {
        //anomalous_reflection Basic transmission loss associated with anomalous propagation
        //   This function computes the basic transmission loss associated with
        //   anomalous propagation as defined in ITU-R P.2001-3 (Attachment D)
        //
        //     Input parameters:
        //     f       -   Frequency GHz
        //     d       -   Vector of distances di of the i-th profile point (km)
        //     z       -   Radio climatic zones 1 = Sea, 3 = Coastal inland, 4 = Inland
        //                 Vectors d and z each contain n+1 profile points
        //     hts     -   Transmitter antenna height in meters above sea level (i=0)
        //     hrs     -   Receiver antenna height in meters above sea level (i=n)
        //     htea    -   Effective height of Tx antenna above smooth surface(m amsl)
        //     hrea    -   Effective height of Rx antenna above smooth surface (m amsl)
        //     hm      -   Path roughness parameter (m)
        //     thetat  -   Tx horizon elevation angle relative to the local horizontal (mrad)
        //     thetar  -   Rx horizon elevation angle relative to the local horizontal (mrad)
        //     dlt     -   Tx to horizon distance (km)
        //     dlr     -   Rx to horizon distance (km)
        //     phimn   -   Mid-point latitude (deg)
        //     omega   -   the fraction of the path over sea
        //     ae      -   Effective Earth radius (km)
        //     p       -   Percentage of average year for which predicted basic loss
        //                 is not exceeded (%)
        //     q       -   Percentage of average year for which predicted basic loss
        //                 is exceeded (%)
        //
        //     Output parameters:
        //     Lba    -   Basic transmission loss associated with ducting (dB)
        //     Aat    -   Time-dependent loss (dB) c.f. Attachment D.7
        //     Aad    -   Angular distance dependent loss (dB) c.f. Attachment D.6
        //     Aac    -   Total coupling loss to the anomalous propagation mechanism (dB) c.f. Attachment D.5
        //     Dct    -   Coast distance from Tx (km)
        //     Dcr    -   Coast distance from Rx (km)
        //     Dtm    -   Longest continuous land (inland or coastal) section of the path (km)
        //     Dlm    -   Longest continuous inland section of the path (km)
        //
        //     Example:
        //     [Lba, Aat, Aad, Aac, Dct, Dcr, Dtm, Dlm] = tl_anomalous_reflection(f, d, h, z, hts, hrs, htea, hrea, thetat, thetar, dlt, dlr, omega, ae, p, q)

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    15JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version
        //     v2    28OCT19     Ivica Stevanovic, OFCOM         Changes in angular distance dependent loss according to ITU-R P.2001-3
        //
        int nl = d.length;
        double dt = d[nl - 1] - d[0];

        // Attachment D: Anomalous/layer-reflection model

        // D.1 Characterize the radio-climatic zones dominating the path

        // Longest continuous land (inland or coastal) section of the path
        int zoner = 34; // (3 = coastal) + (4 = inland)
        double dtm = longest_cont_dist(d, z, zoner);

        // Longest continuous inland section of the path
        zoner = 4;
        double dlm = longest_cont_dist(d, z, zoner);

        // D.2 Point incidence of ducting

        // Calculate a parameter depending on the longest inland section of the path
        // (D.2.1)

        double tau = 1 - Math.exp(-0.000412 * Math.pow(dlm, 2.41));

        // Calculate parameter mu_1 characterizing the degree to which the path is
        // over land (D.2.2)

        double mu1 = Math.pow((Math.pow(10, (-dtm / (16 - 6.6 * tau))) + Math.pow(10, (-(2.48 + 1.77 * tau)))), 0.2);

        if (mu1 > 1) {
            mu1 = 1;
        }

        // Calculate parameter mu4 given by (D.2.3)
        double mu4;

        if (Math.abs(phimn) <= 70) {

            mu4 = Math.pow(10, ((-0.935 + 0.0176 * Math.abs(phimn)) * Math.log10(mu1)));

        } else {

            mu4 = Math.pow(10, (0.3 * Math.log10(mu1)));

        }

        // The point incidence of anomalous propagation for the path centre (D.2.4)
        double b0;
        if (Math.abs(phimn) <= 70) {

            b0 = mu1 * mu4 * Math.pow(10, (-0.015 * Math.abs(phimn) + 1.67));

        } else {

            b0 = 4.17 * mu1 * mu4;

        }


        // D.3 Site-shielding losses with respect to the anomalous propagatoin mechanism

        // Corrections to Tx and Rx horizon elevation angles (D.3.1)

        double gtr = 0.1 * dlt;
        double grr = 0.1 * dlr;

        // Modified transmitter and receiver horizon elevation angles (D.3.2)

        double thetast = thetat - gtr;  //mrad
        double thetasr = thetar - grr;  //mrad

        // Tx and Rx site-shielding losses with respect to the duct (D.3.3)-(D.3.4)
        double Ast;
        if (thetast > 0) {
            Ast = 20 * Math.log10(1 + 0.361 * thetast * Math.sqrt(f * dlt)) + 0.264 * thetast * Math.pow(f, (1.0 / 3.0));
        } else {
            Ast = 0;
        }

        double Asr;
        if (thetasr > 0) {

            Asr = 20 * Math.log10(1 + 0.361 * thetasr * Math.sqrt(f * dlr)) + 0.264 * thetasr * Math.pow(f, (1.0 / 3.0));

        } else {

            Asr = 0;

        }


        // D.4 Over-sea surface duct coupling corrections

        // Obtain the distance from each terminal to the sea in the direction of the
        // other terminal (D.4.1)

        double[] dts = distance_to_sea(d, z);
        double dct = dts[0];
        double dcr = dts[1];

        // The over-sea surface duct coupling corrections for Tx and Rx
        // (D.4.2)-(D.4.3)

        double Act = 0;
        double Acr = 0;

        if (omega >= 0.75) {

            if (dct <= dlt && dct <= 5) {

                Act = -3 * Math.exp(-0.25 * dct * dct) * (1 + Math.tanh(0.07 * (50 - hts)));

            }

            if (dcr <= dlr && dcr <= 5) {

                Acr = -3 * Math.exp(-0.25 * dcr * dcr) * (1 + Math.tanh(0.07 * (50.0 - hrs)));

            }

        }

        // D.5 Total coupling loss to the anomalous propagation mechanism

        // Empirical correction to account for the increasing attenuation with
        // wavelength in ducted propagation (D.5.2)

        double Alf = 0;

        if (f < 0.5) {

            Alf = (45.375 - 137.0 * f + 92.5 * f * f) * omega;

        }

        // Total coupling losses between the antennas and the anomalous propagation
        // mechanism (D.5.1)

        double Aac = 102.45 + 20 * Math.log10(f * (dlt + dlr)) + Alf + Ast + Asr + Act + Acr;

        // D.6 Angular-distance dependent loss

        // Specific angular attenuation (D.6.1)

        double gammad = 5e-5 * ae * Math.pow(f, (1.0 / 3.0));

        // Adjusted Tx and Rx horizon elevation angles (D.6.2)

        double theta_at = Math.min(thetat, gtr);  // mrad

        double theta_ar = Math.min(thetar, grr);  // mrad

        // Adjusted total path angular distance (D.6.3)

        double theta_a = 1000.0 * dt / ae + theta_at + theta_ar; // mrad

        // Angular-distance dependent loss (D.6.4a,b)

        double Aad = 0;
        if (theta_a > 0) {
            Aad = gammad * theta_a;
        }


        // D.7 Distance and time-dependent loss

        // Distance adjusted for terrain roughness factor (D.7.1)

        double dar = Math.min(dt - dlt - dlr, 40.0);

        // Terrain roughness factor (D.7.2)
        double mu3;

        if (hm > 10) {

            mu3 = Math.exp(-4.6e-5 * (hm - 10.0) * (43.0 + 6.0 * dar));

        } else {

            mu3 = 1;

        }

        // A term required for the path geometr ycorrection  (D.7.3)

        double alpha = -0.6 - 3.5e-9 * Math.pow(dt, 3.1) * tau;

        if (alpha < -3.4) {

            alpha = -3.4;

        }

        // Path geometry factor (D.7.4)

        double mu2 = 500 * dt * dt / (ae * Math.pow((Math.sqrt(htea) + Math.sqrt(hrea)), 2.0));
        mu2 = Math.pow(mu2, alpha);

        if (mu2 > 1) {
            mu2 = 1;
        }

        // Time percentage associated with anomalous propagation adjusted for
        // general location and specific properties of the path (D.7.5)

        double bduct = b0 * mu2 * mu3;

        // An exponent required fo rthe time-dependent loss (D.7.6)

        double Gamma = 1.076 * Math.exp(-1e-6 * Math.pow(dt, 1.13) * (9.51 - 4.8 * Math.log10(bduct) + 0.198 * Math.pow((Math.log10(bduct)), 2.0))) /
                (Math.pow((2.0058 - Math.log10(bduct)), 1.012));

        // Time dependent loss (D.7.7)

        double Aat = -12.0 + (1.2 + 0.0037 * dt) * Math.log10(p / bduct) + 12.0 * Math.pow((p / bduct), Gamma) + 50.0 / q;

        // D.8 Basic transmission loss associated with ducting (D.8.1)

        double Lba = Aac + Aad + Aat;

        double[] out = new double[8];
        out[0] = Lba;
        out[1] = Aat;
        out[2] = Aad;
        out[3] = Aac;
        out[4] = dct;
        out[5] = dcr;
        out[6] = dtm;
        out[7] = dlm;
        return out;
    }


    public double longest_cont_dist(double[] d, int[] zone, int zone_r) {
        //longest_cont_dist Longest continuous path belonging to the zone_r
        //     dm = longest_cont_dist(d, zone, zone_r)
        //     This function computes the longest continuous section of the
        //     great-circle path (km) for a given zone_r
        //
        //     Input arguments:
        //     d       -   vector of distances in the path profile
        //     zone    -   vector of zones in the path profile
        //     zone_r  -   reference zone for which the longest continuous section
        //                 is computed
        //
        //     Output arguments:
        //     dm      -   the longest continuous section of the great-circle path (km) for a given zone_r
        //
        //     Example:
        //     dm = longest_cont_dist(d, zone, zone_r)
        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v1    17NOV16     Ivica Stevanovic, OFCOM         First implementation in Java
        //     v2    18JUL17     Ivica Stevanovic, OFCOM         modified mapping to  GlobCover data format
        //                                                       before: 2 - Inland, 1 - Coastal land, 3 - Sea
        //                                                       now:    4 - Inland, 3 - Coastal land, 1 - Sea


        double dm = 0;
        double dmc = 0;
        double delta;
        int n = d.length;
        if (zone_r == 34) { // inland + coastal land
            for (int i = 0; i < n; i++) {
                if ((zone[i] == 3) || (zone[i] == 4)) {
                    if (i == 0) {
                        delta = (d[1] - d[0]) / 2.0;
                    } else if (i == n - 1) {
                        delta = (d[n - 1] - d[n - 2]) / 2.0;
                    } else {
                        delta = (d[i + 1] - d[i - 1]) / 2.0;
                    }

                    dmc = dmc + delta;
                    dm = Math.max(dm, dmc);
                } else {
                    dmc = 0;
                }
            }
            return dm;
        } else {
            for (int i = 0; i < n; i++) {
                if (zone[i] == zone_r) {
                    if (i == 0) {
                        delta = (d[1] - d[0]) / 2.0;
                    } else if (i == n - 1) {
                        delta = (d[n - 1] - d[n - 2]) / 2.0;
                    } else {
                        delta = (d[i + 1] - d[i - 1]) / 2.0;
                    }

                    dmc = dmc + delta;
                    dm = Math.max(dm, dmc);
                } else {
                    dmc = 0;
                }
            }
            return dm;
        }
    }

    public double[] distance_to_sea(double[] d, int[] zone) {
        //distance_to_sea Distance to the sea in the direction of the other terminal
        //     This function computes the distance from eacht terminal to the sea in
        //     the direction of the other terminal
        //
        //     Input arguments:
        //     d       -   vector of distances in the path profile
        //     zone    -   vector of zones in the path profile
        //                 4 - Inland, 3 - Coastal, 1 - See
        //
        //     Output arguments:
        //     dct      -  coast distance from transmitter
        //     dcr      -  coast distance from receiver
        //
        //     Example:
        //     [dct, dcr] = distance_to_sea(d, zone)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    15JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    18JUL17     Ivica Stevanovic, OFCOM         Initial Java Version

        // Find the path profile points belonging to sea (zone = 1) starting from Tx

        int nl = d.length;
        double dct;
        double dcr;
        int nt = -1;
        int nr = -1;

        for (int i = 0; i < nl; i++) {
            if (zone[i] == 1) {
                nt = i;
                break;
            }
        }

        for (int i = nl - 1; i >= 0; i--) {
            if (zone[i] == 1) {
                nr = i;
                break;
            }
        }


        if (nt == -1 && nr == -1) { // complete path is inland or inland coast
            dct = d[nl - 1];
            dcr = d[nl - 1];

        } else {

            if (nt == 0) { // Tx is over sea
                dct = 0;
            } else {
                dct = (d[nt] + d[nt - 1]) / 2 - d[0];
            }

            if (nr == (nl - 1)) { // Rx is over sea
                dcr = 0;
            } else {
                dcr = d[nl - 1] - (d[nr] + d[nr + 1]) / 2;
            }
        }


        double[] out = new double[2];
        out[0] = dct;
        out[1] = dcr;
        return out;
    }

    public double[] tl_troposcatter(P2001DigitalMaps maps, double f, double dt, double thetat, double thetar, double thetae, double phicvn, double phicve, double phitn, double phite, double phirn, double phire, double Gt, double Gr, double ae, double p) {
        //tl_troposcatter Troposcatter basic transmission loss
        //   This function computes the troposcatter basic transmission loss
        //   as defined in ITU-R P.2001-3 (Attachment E)
        //
        //     Input parameters:
        //     maps    -   Object containing Digital Maps
        //     f       -   Frequency GHz
        //     dt      -   Total distance (km)
        //     thetat  -   Tx horizon elevation angle relative to the local horizontal (mrad)
        //     thetar  -   Rx horizon elevation angle relative to the local horizontal (mrad)
        //     thetae  -   Angle subtended by d km at centre of spherical Earth (rad)
        //     phicvn  -   Troposcatter common volume latitude (deg)
        //     phicve  -   Troposcatter common volume longitude (deg)
        //     phitn   -   Tx latitude (deg)
        //     phite   -   Tx longitude (deg)
        //     phirn   -   Rx latitude (deg)
        //     phire   -   Rx longitude (deg)
        //     Gt, Gr  -   Gain of transmitting and receiving antenna in the azimuthal direction
        //                 of the path towards the other antenna and at the elevation angle
        //                 above the local horizontal of the other antenna in the case of a LoS
        //                 path, otherwise of the antenna's radio horizon, for median effective
        //                 Earth radius.
        //     ae      -   Effective Earth radius (km)
        //     p       -   Percentage of average year for which predicted basic loss
        //                 is not exceeded (%)
        //
        //     Output parameters:
        //     Lbs    -   Troposcatter basic transmission loss (dB)
        //     theta  -   Scatter angle (mrad)
        //     climzone-  Climate zone (0,1,2,3,4,5,6)
        //
        //     Example:
        //     [Lbs, theta, climzone] = tl_troposcatter(f, dt, thetat, thetar, thetae, phicvn, phicve, phitn, phite, phirn, phire, Gt, Gr, ae, p)

        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    18JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUN17     Ivica Stevanovic, OFCOM         replaced load function calls to increase computational speed
        //     v2    18JUL17     Ivica Stevanovic, OFCOM         Initial Java Version
        //     v3    19JUL17     Ivica Stevanovic, OFCOM         Introduced Digital Maps

        // Attachment E: Troposcatter

        // E.2 Climatic classification

        // Obtain TropoClim for phicvn, phicve from the data file "TropoClim.txt"

        //TropoClim = load('DigitalMaps/TropoClim.txt');
        // The value at the closest grid point to phicvn, phicve should be taken
        int climzone = maps.GetClimZone(phicve, phicvn);


        // In case the troposcatter common volume lies over the sea the climates at
        // both the transmitter and receiver locations should be determined

        if (climzone == 0) {

            int climzoner = maps.GetClimZone(phire, phirn);


            int climzonet = maps.GetClimZone(phite, phitn);

            // if both terminals have a climate zone corresponding to a land point,
            // the climate zone of the path is given by the smaller value of the
            // transmitter and receiver climate zones

            if (climzoner > 0 && climzonet > 0) {
                climzone = Math.min(climzoner, climzonet);

                // if only one terminal has a climate zone corresponding to a land
                // point, then that climate zone defines the climate zone of the path

            } else if (climzonet > 0) {
                climzone = climzonet;

            } else if (climzoner > 0) {
                climzone = climzoner;
            }

        }


        // From Table E.1 assign meteorological and atmospheric parameters M,
        // gamma and equation
        double M, gamma;
        int eq;
        if (climzone == 1) {
            M = 129.6;
            gamma = 0.33;
            eq = 8;

        } else if (climzone == 2) {
            M = 119.73;
            gamma = 0.27;
            eq = 6;

        } else if (climzone == 3) {
            M = 109.3;
            gamma = 0.32;
            eq = 9;

        } else if (climzone == 4) {
            M = 128.5;
            gamma = 0.27;
            eq = 10;

        } else if (climzone == 5) {
            M = 119.73;
            gamma = 0.27;
            eq = 6;

        } else if (climzone == 6) {
            M = 123.2;
            gamma = 0.27;
            eq = 6;

        } else {// climzone == 0
            M = 116;
            gamma = 0.27;
            eq = 7;
        }


        // E.3 Calculation of tropocscatter basic transmission loss

        // The scatter angle (E.1)

        double theta = 1000 * thetae + thetat + thetar; // mrad

        // The loss term dependent on the common vaolume height

        double H = 0.25e-3 * theta * dt;                                                 // Eq (E.3)

        double htrop = 0.125e-6 * theta * theta * ae;                                          // Eq (E.4)

        double LN = 20 * Math.log10(5 + gamma * H) + 4.34 * gamma * htrop;                      // Eq (E.2)

        // Angular distance of the scatter path based on median effective Earth
        // radius

        double ds = 0.001 * theta * ae;                                                  // Eq (E.5)

        // Calculate Y90 (dB) using one od equations (E.6)-(E.10) as selected from
        // table E.1
        double Y90 = 0;

        if (eq == 6) {

            Y90 = -2.2 - (8.1 - 0.23 * Math.min(f, 4)) * Math.exp(-0.137 * htrop);          // Eq (E.6)

        } else if (eq == 7) {

            Y90 = -9.5 - 3 * Math.exp(-0.137 * htrop);                                 // Eq (E.7)

        } else if (eq == 8) {

            if (ds < 100) {
                Y90 = -8.2;

            } else if (ds >= 1000) {
                Y90 = -3.4;

            } else {
                Y90 = 1.006e-8 * Math.pow(ds, 3) - 2.569e-5 * Math.pow(ds, 2) + 0.02242 * ds - 10.2;    // Eq (E.8)
            }

        } else if (eq == 9) {

            if (ds < 100) {
                Y90 = -10.845;

            } else if (ds >= 465) {
                Y90 = -8.4;

            } else {
                Y90 = -4.5e-7 * Math.pow(ds, 3) + 4.45e-4 * Math.pow(ds, 2) - 0.122 * ds - 2.645;       // Eq (E.9)

            }

        } else if (eq == 10) {

            if (ds < 100) {
                Y90 = -11.5;

            } else if (ds >= 550) {
                Y90 = -4;

            } else {
                Y90 = -8.519e-8 * Math.pow(ds, 3) + 7.444e-5 * Math.pow(ds, 2) - 4.18e-4 * ds - 12.1;   // Eq (E.10)

            }

        }

        // Conversion factor given by (E.11)
        double C;
        if (p >= 50) {

            C = 1.26 * Math.pow((-Math.log10((100 - p) / 50)), 0.63);

        } else {
            C = -1.26 * Math.pow((-Math.log10(p / 50)), 0.63);

        }

        // Parameter Yp not exceeded for p% time (E.12)

        double Yp = C * Y90;

        // Limit the value of theta such that theta >= 1e-6

        if (theta < 1e-6) {
            theta = 1e-6;
        }

        // Distance and frequency dependent losses (E.13) and (E.14)

        double Ldist = Math.max(10 * Math.log10(dt) + 30 * Math.log10(theta) + LN, 20 * Math.log10(dt) + 0.573 * theta + 20);

        double Lfreq = 25 * Math.log10(f) - 2.5 * Math.pow((Math.log10(0.5 * f)), 2);

        // Aperture-to-medium copuling loss (E.15)

        double Lcoup = 0.07 * Math.exp(0.055 * (Gt + Gr));

        // Troposcatter basic transmission loss not exceeded for p// time (E.16)

        double Lbs = M + Lfreq + Ldist + Lcoup - Yp;

        double[] out = new double[3];
        out[0] = Lbs;
        out[1] = theta;
        out[2] = climzone;

        return out;
    }

    public double[] gaseous_abs_tropo_t2cv(double rho_sur, double h_sur, double theta_el, double dcv, double f) {
        //gaseous_abs_tropo_t2cv Gaseous absorbtion for tropospheric terminal-common-volume path
        //   This function computes gaseous absorbtion for termina/common-volume
        //   troposcatter path as defined in ITU-R P.2001-3 Attachment F.4
        //
        //   Input parameters:
        //   rho_sur   -   Surface water-vapour density under non rain conditions (g/m^3)
        //   h_sur     -   terrain height (masl)
        //   theta_el  -   elevation angle of path (mrad)
        //   dcv       -   Horizontal distance to the comon volume (km)
        //   f         -   Frequency (GHz), not greater than 54 GHz
        //
        //   Output parameters:
        //   Ao       -   Attenuation due to oxygen (dB)
        //   Aw       -   Attenuation due to water-vapour under non-rain conditions (dB)
        //   Awr      -   Attenuation due to water-vapour under rain conditions (dB)
        //   gamma_o  -   Specific attenuation due to oxygen (dB/km)
        //   gamma_w  -   Specific attenuation due to water-vapour non-rain conditions (dB/km)
        //   gamma_wr -   Specific attenuation due to water-vapour rain conditions (dB/km)
        //
        //   Example:
        //   [Ao, Aw, Awr, gamma_o, gamma_w, gamma_wr] = gaseous_abs_tropo_t2cv(rho_sur, h_sur, theta_el, dcv, f)

        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    18JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUN17     Ivica Stevanovic, OFCOM         Octave compatibility (do -> d0)
        //     v2    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version

        // Use equation (F.6.2) to calculate the sea-level specific attenuation due
        // to water vapour under non-rain conditions gamma_w

        double[] ssla = specific_sea_level_attenuation(f, rho_sur, h_sur);
        double gamma_o = ssla[0];
        double gamma_w = ssla[1];

        // Use equation (F.5.1) to calculate the surface water-vapour density under
        // rain conditions rho_surr

        double rho_surr = water_vapour_density_rain(rho_sur, h_sur);

        // Use equation (F.6.2) to calculate the sea-level specifica ttenuation due
        // to water vapour undr rain conditions gamma_wr

        double[] ssla1 = specific_sea_level_attenuation(f, rho_surr, h_sur);
        double gamma_wr = ssla1[1];

        // Calculate the quantities do and dw for oxygen and water vapour (F.4.1)

        double d0 = 0.65 * Math.sin(0.001 * theta_el) + 0.35 * Math.sqrt(Math.pow((Math.sin(0.001 * theta_el)), 2) + 0.00304);
        d0 = 5 / d0;

        double dw = 0.65 * Math.sin(0.001 * theta_el) + 0.35 * Math.sqrt(Math.pow((Math.sin(0.001 * theta_el)), 2) + 0.00122);
        dw = 2 / dw;

        // Effective distances for oxygen and water vapour (F.4.2)

        double deo = d0 * (1 - Math.exp(-dcv / d0)) * Math.exp(-h_sur / 5000);


        double dew = dw * (1 - Math.exp(-dcv / dw)) * Math.exp(-h_sur / 2000);

        // Attenuations due to oxygen, and for water vapour for both non-rain and
        // rain conditions (F.4.3)

        double Ao = gamma_o * deo;

        double Aw = gamma_w * dew;

        double Awr = gamma_wr * dew;

        double[] out = new double[6];

        out[0] = Ao;
        out[1] = Aw;
        out[2] = Awr;
        out[3] = gamma_o;
        out[4] = gamma_w;
        out[5] = gamma_wr;


        return out;
    }

    public double[] gaseous_abs_tropo(P2001DigitalMaps maps, double phi_te, double phi_tn, double phi_re, double phi_rn, double h1, double hn, double thetatpos, double thetarpos, double dtcv, double drcv, double f) {
        //gaseous_abs_tropo Gaseous absorbtion for a troposcatter path
        //   This function computes gaseous absorbtion for a complete troposcater
        //   path, from Tx to Rx via the common scattering volume
        //   and water-vapour as defined in ITU-R P.2001-3 Attachment F.3
        //   The formulas are valid for frequencies not greater than 54 GHz.
        //
        //   Input parameters:
        //   maps      -   Object containing digital maps
        //   phi_te    -   Tx Longitude (deg)
        //   phi_tn    -   Tx Latitude  (deg)
        //   phi_re    -   Rx Longitude (deg)
        //   phi_rn    -   Rx Latitude  (deg)
        //   h1        -   Ground height at the transmitter point of the profile (masl)
        //   hn        -   Ground height at the receiver point of the profile (masl)
        //   thetatpos -   Horizon elevation angle relative to the local horizontal
        //                 as viewed from Tx (limited to be positive) (mrad)
        //   thetarpos -   Horizon elevation angle relative to the local horizontal
        //                 as viewed from Rx (limited to be positive) (mrad)
        //   dtcv      -   Tx terminal to troposcatter common volume distance
        //   drcv      -   Rx terminal to troposcatter common volume distance
        //   f         -   Frequency (GHz), not greater than 54 GHz
        //
        //   Output parameters:
        //   Aos       -   Attenuation due to oxygen for the complete troposcatter path (dB)
        //   Aws       -   Attenuation due to water-vapour under non-rain conditions for the complete path(dB)
        //   Awrs      -   Attenuation due to water-vapour under rain conditions for the complete path(dB)
        //   Aotcv     -   Attenuation due to oxygen for the Tx-cv path (dB)
        //   Awtcv     -   Attenuation due to water-vapour under non-rain conditions for the Tx-cv path(dB)
        //   Awrtcv    -   Attenuation due to water-vapour under rain conditions for the Tx-cv path(dB)
        //   Aorcv     -   Attenuation due to oxygen for the Rx-cv path (dB)
        //   Awrcv     -   Attenuation due to water-vapour under non-rain conditions for the Rx-cv path(dB)
        //   Awrrcv    -   Attenuation due to water-vapour under rain conditions for the Rx-cv path(dB)
        //   Wvsurtx   -   Surface water-vapour density under non-rain conditions at the Tx (g/m^3)
        //   Wvsurrx   -   Surface water-vapour density under non-rain conditions at the Rx (g/m^3)


        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    18JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUN17     Ivica Stevanovic, OFCOM         replaced load function calls to increase computational speed
        //     v2    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version
        //     v3    19JUL17     Ivica Stevanovic, OFCOM         Introduced Digital Maps

        // Obtain surface water-vapour density under non-rain conditions at the
        // location of Tx from the data file surfwv_50_fixed.txt

        //surfwv_50_fixed = load('DigitalMaps/surfwv_50_fixed.txt');

        // Find rho_sur from file surfwv_50_fixed.txt as a bilinear interpolation

        double rho_sur = maps.GetRhoSur(phi_te, phi_tn);

        double Wvsurtx = rho_sur;


        // Use the method in Attachment F.4 to get the gaseous attenuations due to
        // oxygen and for water vapour under both non-rain and rain conditions for
        // the Tx-cv path (F.3.1)
        double h_sur = h1;
        double theta_el = thetatpos;
        double dcv = dtcv;

        double[] gat2tcv = gaseous_abs_tropo_t2cv(rho_sur, h_sur, theta_el, dcv, f);
        double Aotcv = gat2tcv[0];
        double Awtcv = gat2tcv[1];
        double Awrtcv = gat2tcv[2];

        // Obtain surface water-vapour density under non-rain conditions at the
        // location of Rx from the data file surfwv_50_fixed.txt

        // Find rho_sur from file surfwv_50_fixed.txt as a bilinear interpolation

        rho_sur = maps.GetRhoSur(phi_re, phi_rn);

        double Wvsurrx = rho_sur;


        // Use the method in Attachment F.4 to get the gaseous attenuations due to
        // oxygen and for water vapour under both non-rain and rain conditions for
        // the Rx-cv path (F.3.2)
        h_sur = hn;
        theta_el = thetarpos;
        dcv = drcv;

        double[] gat2rcv = gaseous_abs_tropo_t2cv(rho_sur, h_sur, theta_el, dcv, f);
        double Aorcv = gat2rcv[0];
        double Awrcv = gat2rcv[1];
        double Awrrcv = gat2rcv[2];

        // Gaseous attenuations for the complete troposcatter path (F.3.3)

        double Aos = Aotcv + Aorcv;

        double Aws = Awtcv + Awrcv;

        double Awrs = Awrtcv + Awrrcv;

        double[] out = new double[11];
        out[0] = Aos;
        out[1] = Aws;
        out[2] = Awrs;
        out[3] = Aotcv;
        out[4] = Awtcv;
        out[5] = Awrtcv;
        out[6] = Aorcv;
        out[7] = Awrcv;
        out[8] = Awrrcv;
        out[9] = Wvsurtx;
        out[10] = Wvsurrx;

        return out;
    }


    double[] tl_sporadic_e(P2001DigitalMaps maps, double f, double dt, double thetat, double thetar, double phimn, double phime, double phitn, double phite, double phirn, double phire, double dlt, double dlr, double ae, double Re, double p) {
        //tl_sporadic_e Sporadic-E transmission loss
        //   This function computes the basic transmission loss due to sporadic-E
        //   propagation as defined in ITU-R P.2001-3 (Attachment G)
        //
        //     Input parameters:
        //     maps    -   Object containing digital maps
        //     f       -   Frequency GHz
        //     dt      -   Total distance (km)
        //     thetat  -   Tx horizon elevation angle relative to the local horizontal (mrad)
        //     thetar  -   Rx horizon elevation angle relative to the local horizontal (mrad)
        //     phimn   -   Mid-point latitude (deg)
        //     phime   -   Mid-point longitude (deg)
        //     phitn   -   Tx latitude (deg)
        //     phite   -   Tx longitude (deg)
        //     phirn   -   Rx latitude (deg)
        //     phire   -   Rx longitude (deg)
        //     dlt     -   Tx to horizon distance (km)
        //     dlr     -   Rx to horizon distance (km)
        //     ae      -   Effective Earth radius (km)
        //     Re      -   Average Earth radius (km)
        //     p       -   Percentage of average year for which predicted basic loss
        //                 is not exceeded (%)
        //
        //     Output parameters:
        //     Lbe    -   Basic transmission loss due to sporadic-E propagation (dB)
        //     Lbes1  -   Sporadic-E 1-hop basic transmission loss (dB)
        //     Lbes2  -   Sporadic-E 2-hop basic transmission loss (dB)
        //     Lp1t/r -   Diffraction losses at the two terminals for 1-hop propagation (dB)
        //     Lp2t/r -   Diffraction losses at the two terminals for 2-hop propagation (dB)
        //     Gamma1/2   Ionospheric loss for 1/2 hops (dB)
        //     FoEs1/2hop FoEs for 1/2 hop(s) sporadic-E propagation
        //     Phi1qe -   Longitude of the one-quarter point
        //     Phi1qn -   Latitude of the one-quarter point
        //     Phi3qe -   Longitude of the three-quarter point
        //     Phi3qn -   Latitude of the three-quarter point
        //


        //
        //     Rev   Date        Author                          Description
        //     -------------------------------------------------------------------------------
        //     v0    18JUL16     Ivica Stevanovic, OFCOM         Initial version
        //     v1    13JUN17     Ivica Stevanovic, OFCOM         replaced load function calls to increase computational speed
        //     v2    18JUL17     Ivica Stevanovic, OFCOM         Initial Java version
        //     v3    19JUL17     Ivica Stevnaovic, OFCOM         Introduced Digital Maps

        // Attachment G: Sporadic-E propagation

        // G.2 Derivation of FoEs
        double foes1;
        double foes2;
        double p1;
        double p2;

        if (p < 1) {
            p1 = 0.1;
            p2 = 1;

        } else if (p > 10) {
            p1 = 10;
            p2 = 50;

        } else {
            p1 = 1;
            p2 = 10;

        }

        foes1 = maps.GetFoEs(phime, phimn, p1);
        foes2 = maps.GetFoEs(phime, phimn, p2);


        // G.2 1-hop propagation


        // Find foes1/2 from the correspoinding files for the path mid-point - as a bilinear interpolation

        double FoEs1hop = foes1 + (foes2 - foes1) * (Math.log10(p / p1)) / (Math.log10(p2 / p1));          // Eq (G.1.1)

        // Ionospheric loss for one hop

        double Gamma1 = (40 / (1 + dt / 130 + Math.pow((dt / 250), 2)) + 0.2 * Math.pow((dt / 2600), 2)) * Math.pow((1000 * f / FoEs1hop), 2) + Math.exp((dt - 1660) / 280);     // Eq (G.2.1)

        // Slope path length

        double hes = 120;

        double l1 = 2 * Math.pow((Math.pow(ae, 2) + Math.pow((ae + hes), 2) - 2 * ae * (ae + hes) * Math.cos(dt / (2 * ae))), 0.5); // Eq (G.2.2)

        // free space loss for the slope distance:

        double Lbfs1 = tl_free_space(f, l1);                                             // Eq (G.2.3)

        // ray take-off angle above the horizontal at both terminals for 1 hop

        double alpha1 = dt / (2 * ae);

        double epsr1 = 0.5 * Math.PI - Math.atan(ae * Math.sin(alpha1) / (hes + ae * (1 - Math.cos(alpha1)))) - alpha1;    // Eq (G.2.4)

        // Diffraction angles for the two terminals

        double delta1t = 0.001 * thetat - epsr1;
        double delta1r = 0.001 * thetar - epsr1;                                         // (G.2.5)

        // Diffraction parameters (G.2.6)
        double nu1t, nu1r;

        if (delta1t >= 0) {

            nu1t = 3.651 * Math.sqrt(1000 * f * dlt * (1 - Math.cos(delta1t)) / Math.cos(0.001 * thetat));

        } else {

            nu1t = -3.651 * Math.sqrt(1000 * f * dlt * (1 - Math.cos(delta1t)) / Math.cos(0.001 * thetat));

        }


        if (delta1r >= 0) {

            nu1r = 3.651 * Math.sqrt(1000 * f * dlr * (1 - Math.cos(delta1r)) / Math.cos(0.001 * thetar));

        } else {

            nu1r = -3.651 * Math.sqrt(1000 * f * dlr * (1 - Math.cos(delta1r)) / Math.cos(0.001 * thetar));

        }

        // Diffraction lossess at the two terminals (G.2.7)

        double Lp1t = dl_knife_edge(nu1t);

        double Lp1r = dl_knife_edge(nu1r);

        // Sporadic-E 1-hop basic transmission loss (G.2.8)

        double Lbes1 = Lbfs1 + Gamma1 + Lp1t + Lp1r;

        // G.3 2-hop propagation

        // Latitude and longitude of the one-quarter point

        double d1q = 0.25 * dt;

        double[] gcp = great_circle_path(phire, phite, phirn, phitn, Re, d1q);

        double Phi1qe = gcp[0];
        double Phi1qn = gcp[1];

        // Latitude and longitude of the one-quarter point

        double d3q = 0.75 * dt;

        double[] gcp3 = great_circle_path(phire, phite, phirn, phitn, Re, d3q);

        double Phi3qe = gcp3[0];
        double Phi3qn = gcp3[1];

        // foes for one-quarter point
        // Map phime (-180, 180) to loncnt (0,360);

        double phie = Phi1qe;
        double phin = Phi1qn;

        foes1 = maps.GetFoEs(phie, phin, p1);
        foes2 = maps.GetFoEs(phie, phin, p2);

        double FoEs2hop1q = foes1 + (foes2 - foes1) * (Math.log10(p / p1)) / (Math.log10(p2 / p1));          // Eq (G.1.1)

        // foes for three-quarter point
        // Map phie (-180, 180) to loncnt (0,360);

        phie = Phi3qe;
        phin = Phi3qn;

        foes1 = maps.GetFoEs(phie, phin, p1);
        foes2 = maps.GetFoEs(phie, phin, p2);

        double FoEs2hop3q = foes1 + (foes2 - foes1) * (Math.log10(p / p1)) / (Math.log10(p2 / p1));          // Eq (G.1.1)

        // Obtain FoEs2hop as the lower of the two values calculated above

        double FoEs2hop = Math.min(FoEs2hop1q, FoEs2hop3q);

        // Ionospheric laps for two hops (G.3.1)

        double Gamma2 = (40 / (1 + (dt / 260) + Math.pow((dt / 500), 2)) + 0.2 * Math.pow((dt / 5200), 2)) * Math.pow((1000 * f / FoEs2hop), 2) + Math.exp((dt - 3220) / 560);

        // Slope path length

        double l2 = 4 * Math.pow((ae * ae + Math.pow((ae + hes), 2) - 2 * ae * (ae + hes) * Math.cos(dt / (4 * ae))), 0.5);       // Eq (G.3.2)

        // Free-space loss for this slope

        double Lbfs2 = tl_free_space(f, l2);                                             // Eq (G.3.3)

        // Ray take-off angle above the local horiozntal at both terminals for 2
        // hops (G.3.4)

        double alpha2 = dt / (4 * ae);

        double epsr2 = 0.5 * Math.PI - Math.atan(ae * Math.sin(alpha2) / (hes + ae * (1 - Math.cos(alpha2)))) - alpha2;

        // Diffraction angles for the two terminals (G.3.5)

        double delta2t = 0.001 * thetat - epsr2;
        double delta2r = 0.001 * thetar - epsr2;

        // Corresponding diffraction parameters (G.3.6)
        double nu2t, nu2r;

        if (delta2t >= 0) {

            nu2t = 3.651 * Math.sqrt(1000 * f * dlt * (1 - Math.cos(delta2t)) / Math.cos(0.001 * thetat));

        } else {

            nu2t = -3.651 * Math.sqrt(1000 * f * dlt * (1 - Math.cos(delta2t)) / Math.cos(0.001 * thetat));

        }


        if (delta2r >= 0) {

            nu2r = 3.651 * Math.sqrt(1000 * f * dlr * (1 - Math.cos(delta2r)) / Math.cos(0.001 * thetar));

        } else {

            nu2r = -3.651 * Math.sqrt(1000 * f * dlr * (1 - Math.cos(delta2r)) / Math.cos(0.001 * thetar));

        }

        // Diffraction lossess at the two terminals (G.3.7)

        double Lp2t = dl_knife_edge(nu2t);

        double Lp2r = dl_knife_edge(nu2r);

        // Sporadic-E two-hop basic transmission loss

        double Lbes2 = Lbfs2 + Gamma2 + Lp2t + Lp2r;

        // G.4 Basic transmission loss (G.4.1)
        double Lbe;

        if (Lbes1 < Lbes2 - 20) {
            Lbe = Lbes1;
        } else if (Lbes2 < Lbes1 - 20) {
            Lbe = Lbes2;
        } else {
            Lbe = -10 * Math.log10(Math.pow(10, (-0.1 * Lbes1)) + Math.pow(10, (-0.1 * Lbes2)));
        }

        double[] out = new double[15];

        out[0] = Lbe;
        out[1] = Lbes1;
        out[2] = Lbes2;
        out[3] = Lp1t;
        out[4] = Lp2t;
        out[5] = Lp1r;
        out[6] = Lp2r;
        out[7] = Gamma1;
        out[8] = Gamma2;
        out[9] = FoEs1hop;
        out[10] = FoEs2hop;
        out[11] = Phi1qe;
        out[12] = Phi1qn;
        out[13] = Phi3qe;
        out[14] = Phi3qn;

        return out;
    }

}
