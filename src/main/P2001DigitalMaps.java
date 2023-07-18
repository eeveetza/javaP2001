package main;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


public class P2001DigitalMaps {


    private DigitalMap _mapSdN;
    private DigitalMap _mapSdNsup;
    private DigitalMap _mapSdNsub;
    private DigitalMap _mapDndz_01;
    private DigitalMap _mapRho_sur;
    private DigitalMap _mapEsarain_Pr6;
    private DigitalMap _mapEsarain_Mt;
    private DigitalMap _mapEsarain_Beta;
    private DigitalMap _mapData_h0;
    private DigitalMap _mapFoEs50;
    private DigitalMap _mapFoEs10;
    private DigitalMap _mapFoEs01;
    private DigitalMap _mapFoEs0_1;
    private DigitalMap _mapTropoClim;
    private DigitalMap _mapDN50;
    private DigitalMap _mapN050;


    public P2001DigitalMaps() {

        _mapSdN = new DigitalMap("DN_Median.txt", true);
        _mapSdNsup = new DigitalMap("DN_SupSlope.txt", true);
        _mapSdNsub = new DigitalMap("DN_SubSlope.txt", true);
        _mapDndz_01 = new DigitalMap("dndz_01.txt", true);
        _mapRho_sur = new DigitalMap("surfwv_50_fixed.txt", true);
        _mapEsarain_Pr6 = new DigitalMap("Esarain_Pr6_v5.txt", true);
        _mapEsarain_Mt = new DigitalMap("Esarain_Mt_v5.txt", true);
        _mapEsarain_Beta = new DigitalMap("Esarain_Beta_v5.txt", true);
        _mapData_h0 = new DigitalMap("h0.txt", true);
        _mapFoEs50 = new DigitalMap("FoEs50.txt", true);
        _mapFoEs10 = new DigitalMap("FoEs10.txt", true);
        _mapFoEs01 = new DigitalMap("FoEs01.txt", true);
        _mapFoEs0_1 = new DigitalMap("FoEs0.1.txt", true);
        _mapDN50 = new DigitalMap("DN50.TXT", true);
        _mapN050 = new DigitalMap("N050.TXT", true);
        _mapTropoClim = new DigitalMap("TropoClim.txt", false);

    }

    public double GetSdN(double lon, double lat) {
        return _mapSdN.GetInterpolatedValue(lon, lat);
    }

    public double GetSdNsup(double lon, double lat) {
        return _mapSdNsup.GetInterpolatedValue(lon, lat);
    }

    public double GetSdNsub(double lon, double lat) {
        return _mapSdNsub.GetInterpolatedValue(lon, lat);
    }

    public double GetNd65m1(double lon, double lat) {
        return _mapDndz_01.GetInterpolatedValue(lon, lat);
    }

    public double GetRhoSur(double lon, double lat) {
        return _mapRho_sur.GetInterpolatedValue(lon, lat);
    }

    public double GetPr6(double lon, double lat) {
        return _mapEsarain_Pr6.GetInterpolatedValue(lon, lat);
    }

    public double GetMt(double lon, double lat) {
        return _mapEsarain_Mt.GetInterpolatedValue(lon, lat);
    }

    public double GetBetaRain(double lon, double lat) {
        return _mapEsarain_Beta.GetInterpolatedValue(lon, lat);
    }

    public double GetH0(double lon, double lat) {
        return _mapData_h0.GetInterpolatedValue(lon, lat);
    }

    public double GetDN50(double lon, double lat) { return _mapDN50.GetInterpolatedValue(lon, lat);  }

    public double GetN050(double lon, double lat) { return _mapN050.GetInterpolatedValue(lon, lat);  }

    public double GetFoEs(double lon, double lat, double p) {
        double foEs;
        if (p == 50) {
            foEs = _mapFoEs50.GetInterpolatedValue(lon, lat);
        } else if (p == 10) {
            foEs = _mapFoEs10.GetInterpolatedValue(lon, lat);
        } else if (p == 1) {
            foEs = _mapFoEs01.GetInterpolatedValue(lon, lat);
        } else if (p == 0.1) {
            foEs = _mapFoEs0_1.GetInterpolatedValue(lon, lat);
        } else {
            throw new RuntimeException("Allowed values for p: 0.1, 1, 10, 50.");
        }
        return foEs;
    }

    public int GetClimZone(double lon, double lat) {
        return (int) _mapTropoClim.GetClosestGridPointValue(lon, lat);
    }


    class DigitalMap {
        double[][] _map;
        int _sizeX, _sizeY;
        private double _spacing;

        public DigitalMap(String path, boolean lastAndFirstColumnEqual) {

            String name = "/maps/" + path;
            try {
                InputStream inputStream = getClass().getResourceAsStream(name);

                List<String> lines = new ArrayList<>();

                BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
                String line;
                while (null != (line = br.readLine())) {
                    lines.add(line);
                }

                _sizeY = lines.size();

                _sizeX = ParseLine(lines.get(0)).length;

                _map = new double[_sizeY][_sizeX];

                for (int i = 0; i < _sizeY; i++) { /* DO */
                    double[] data = ParseLine(lines.get(i));
                    for (int j = 0; j < _sizeX; j++) {
                        _map[i][j] = data[j];
                    }

                }

            } catch (Exception ex) {
                throw new IllegalArgumentException("Could not load map: '" + name + "'");
            }

            if (lastAndFirstColumnEqual) {
                _spacing = 360.0 / (_map[1].length - 1);
            } else {
                _spacing = 360.0 / (_map[1].length);
            }
        }

        private double[] ParseLine(String line) {
            String[] parts = line.trim().split("\\s+");
            ;

            double[] data = new double[parts.length];

            for (int i = 0; i < parts.length; i++) {
                data[i] = Double.parseDouble(parts[i]);
                //System.out.printf(     "parts.length    =  %f\n"  ,data[i]);
            }
            return data;
        }

        public double GetClosestGridPointValue(double longitude, double latitude) {
            double longitudeOffset = longitude + 180.0;
            double latitudeOffset = 90.0 - latitude;
            int latitudeIndex = (int) (latitudeOffset / _spacing);
            int longitudeIndex = (int) (longitudeOffset / _spacing);

            latitudeIndex %= _sizeY;
            longitudeIndex %= _sizeX;


            double val = _map[latitudeIndex][longitudeIndex];
            return val;
        }

        public double GetInterpolatedValue(double longitude, double latitude) {
            double longitudeOffset = longitude;

            if (longitude < 0.0) {
                longitudeOffset = longitude + 360.0;
            }

            double latitudeOffset = 90.0 - latitude;
            int latitudeIndex = (int) (latitudeOffset / _spacing);
            int longitudeIndex = (int) (longitudeOffset / _spacing);
            //System.out.printf(     "latitudeIndex    =  %d\n"  ,latitudeIndex);
            //System.out.printf(     "longitudeIndex    =  %d\n"  ,longitudeIndex);
            double latitudeFraction = (latitudeOffset / _spacing) - latitudeIndex;
            double longitudeFraction = (longitudeOffset / _spacing) - longitudeIndex;

            double value_ul = _map[latitudeIndex][longitudeIndex];
            double value_ur = _map[latitudeIndex][(longitudeIndex + 1) % _sizeX];
            double value_ll = _map[(latitudeIndex + 1) % _sizeY][longitudeIndex];
            double value_lr = _map[(latitudeIndex + 1) % _sizeY][(longitudeIndex + 1) % _sizeX];

            double interpolatedHeight1 = (longitudeFraction * (value_ur - value_ul)) + value_ul;
            double interpolatedHeight2 = (longitudeFraction * (value_lr - value_ll)) + value_ll;
            double interpolatedHeight3 = latitudeFraction * (interpolatedHeight2 - interpolatedHeight1) + interpolatedHeight1;

            return interpolatedHeight3;
        }
    }
}

