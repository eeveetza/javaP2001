package test;

import main.P2001;
import main.P2001DigitalMaps;

import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;


public class P2001Test {
    // the results are compared to the MATLAB implementation of Recommendation ITU-R P.2001-6
    // the test is passed when the results for transmission loss are within 0.02 dB of difference
    // different antenna heights, time percentages, frequencies, distances, and profiles are tested

    //     Rev   Date        Author                          Description
    //     -------------------------------------------------------------------------------
    //     v0    19JUL17     Ivica Stevanovic, OFCOM         Initial version in Java
    //     v1    29OCT19     Ivica Stevanovic, OFCOM         Changes in angular distance dependent loss according to ITU-R P.2001-3
    //     v2    13JUL21     Ivica Stevanovic, OFCOM         Changes in free-space loss according to ITU-R P.2001-4
    //     v3    15NOV25     Ivica Stevanovic, OFCOM         Updated version number only (P.2001-6, the source identical to P.2001-4
    //                                                       Introdiced reading and handing profile/results in .csv

    TestUtil util;

    @Before
    public void setup() {

        util = new TestUtil(0.02);

    }


    @Test
    public void test1() {
        // test for Digital Maps
        P2001DigitalMaps maps = new P2001DigitalMaps();

        // GetSdN:

        double lon = 0;
        double lat = 0;
        double sdn = maps.GetSdN(lon, lat);
        util.assertDoubleEquals(53.04800, sdn);

        lon = 360;
        lat = 90;
        sdn = maps.GetSdN(lon, lat);
        util.assertDoubleEquals(40.725999999999999, sdn);

        lon = 360;
        lat = -90;
        sdn = maps.GetSdN(lon, lat);
        util.assertDoubleEquals(38.450000000000003, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetSdN(lon, lat);
        util.assertDoubleEquals(48.875640189157593, sdn);

        //GetSdNsup

        lon = 0;
        lat = 0;
        sdn = maps.GetSdNsup(lon, lat);
        util.assertDoubleEquals(23.020000000000000, sdn);

        lon = 360;
        lat = 90;
        sdn = maps.GetSdNsup(lon, lat);
        util.assertDoubleEquals(7.781000000000000, sdn);

        lon = 360;
        lat = -90;
        sdn = maps.GetSdNsup(lon, lat);
        util.assertDoubleEquals(9.683999999999999, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetSdNsup(lon, lat);
        util.assertDoubleEquals(18.607236359102522, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetSdNsub(lon, lat);
        util.assertDoubleEquals(13.873727681679060, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetNd65m1(lon, lat);
        util.assertDoubleEquals(-244.7147699169856, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetRhoSur(lon, lat);
        util.assertDoubleEquals(2.566024595810121, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetPr6(lon, lat);
        util.assertDoubleEquals(16.148936372289601, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetMt(lon, lat);
        util.assertDoubleEquals(332.6332955435777, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetBetaRain(lon, lat);
        util.assertDoubleEquals(0.151614312501983, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetH0(lon, lat);
        util.assertDoubleEquals(3.672142137728531, sdn);


        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetFoEs(lon, lat, 50);
        util.assertDoubleEquals(3.332369311947466, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetFoEs(lon, lat, 10);
        util.assertDoubleEquals(5.282885057803198, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetFoEs(lon, lat, 1);
        util.assertDoubleEquals(9.496184655973734, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetFoEs(lon, lat, 0.1);
        util.assertDoubleEquals(13.975835529041827, sdn);

        lon = -69.4801969319433;
        lat = -36.0460523359984;
        sdn = maps.GetClimZone(lon, lat);
        util.assertDoubleEquals(4, sdn);

        lon = 53.222;
        lat = 87.046;
        sdn = maps.GetClimZone(lon, lat);
        util.assertDoubleEquals(0, sdn);

    }


    @Test
    public void test2() {
        // Reference validation examples from ITU-R 
        P2001DigitalMaps maps = new P2001DigitalMaps();

        P2001 p2001calculator = new P2001();

        int sizeY = 0;

        // path to the directory where profiles are located

        String directoryPath = "src/test/validation_examples/profiles/";

        // Using File class create an object for specific directory
        File directory = new File(directoryPath);

        // Using listFiles method we get all the files of a directory
        // return type of listFiles is array
        File[] files = directory.listFiles();

        // Get name of the all files present in that path
        if (files != null) {
            for (File file : files) {
                //System.out.println(file.getName());

                List<String> lines = new ArrayList<>();

                String file_rel = directoryPath + file.getName();

                // read all the lines from the profile
                try {

                    InputStream inputStream = new FileInputStream(file_rel);
                    InputStreamReader inputStreamReader = new InputStreamReader(inputStream);
                    BufferedReader br = new BufferedReader(inputStreamReader);
                    String line;
                    while (null != (line = br.readLine())) {

                        lines.add(line);

                    }

                    sizeY = lines.size();
                    inputStream.close();

                    double[] d = new double[sizeY - 9];
                    double[] h = new double[sizeY - 9];
                    int[] zone = new int[sizeY - 9];

                    // extract the vectors for distance, terrain height, clutter height and zone
                    // noting that the first line is the table header

                    for (int i = 9; i < sizeY; i++) { /* DO */

                        String[] parts = lines.get(i).trim().split(",");
                        d[i - 9] = Double.parseDouble(parts[0]);
                        h[i - 9] = Double.parseDouble(parts[1]);
                        zone[i - 9] = Integer.parseInt(parts[2]);

                    }

                    // at this point, the profile is read
                    // read the reference results next

                    List<String> lines_r = new ArrayList<>();
                    String directoryPath_r = "src/test/validation_examples/results/";

                    String file_rel_r = directoryPath_r + file.getName().replaceFirst("profile", "results");
                    //System.out.println(file_rel_r);

                    InputStream inputStream_r = new FileInputStream(file_rel_r);
                    InputStreamReader inputStreamReader_r = new InputStreamReader(inputStream_r);
                    BufferedReader br_r = new BufferedReader(inputStreamReader_r);
                    String line_r;
                    while (null != (line_r = br_r.readLine())) {

                        lines_r.add(line_r);
                    }

                    sizeY = lines_r.size();
                    inputStream_r.close();

                    for (int i = 1; i < sizeY; i++) { /* DO */

                        String[] parts = lines_r.get(i).trim().split(",");

                        int FlagVP = Integer.parseInt(parts[0]);
                        double GHz = Double.parseDouble(parts[1]);

                        double Grx = Double.parseDouble(parts[2]);
                        double Gtx = Double.parseDouble(parts[3]);
                        double Hrg = Double.parseDouble(parts[4]);
                        double Htg = Double.parseDouble(parts[5]);
                        double Phire = Double.parseDouble(parts[6]);
                        double Phirn = Double.parseDouble(parts[7]);

                        double Phite = Double.parseDouble(parts[8]);
                        double Phitn = Double.parseDouble(parts[9]);
                        double Tpc = Double.parseDouble(parts[10]);

                        double Lb_ref = Double.parseDouble(parts[78]);

                        double Lb = p2001calculator.tl_p2001(maps, d, h, zone, GHz, Tpc, Phire, Phirn, Phite, Phitn, Hrg, Htg, Grx, Gtx, FlagVP);
                        util.assertDoubleEquals(Lb_ref, Lb);

                    }


                } catch (Exception ex) {

                    throw new IllegalArgumentException("Could not load the file: '" + file_rel + "'");
                }


            }


        } else {
            System.out.println("Did not find any files");
        }


    }

}




