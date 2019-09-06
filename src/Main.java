import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.similarity.Tanimoto;
import java.util.concurrent.TimeUnit;

import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.*;
import java.nio.charset.StandardCharsets;

public class Main {

    public static int nombre = 10000;

    public static void main(String[] args) {



        try{

            //input file
            InputStream flux=new FileInputStream(new File("molecules.data"));
            InputStreamReader lecture=new InputStreamReader(flux);
            BufferedReader buff=new BufferedReader(lecture);

            //output file 1- similarity
            File fic = new File("tanimoto_fcfp6.result");
            FileOutputStream fichier = new FileOutputStream(fic);
            BufferedWriter res = new BufferedWriter(new OutputStreamWriter(fichier));


            //output file 2- fingerprint computation time

            File fich1 = new File("empreintes_fcfp6.temps");
            FileOutputStream f1 = new FileOutputStream(fich1);
            BufferedWriter res1 = new BufferedWriter(new OutputStreamWriter(f1));


            // output file 3 - similarity computation time

            File fich2 = new File("tanimoto_fcfp6.temps");
            FileOutputStream f2 = new FileOutputStream(fich2);
            BufferedWriter res2 = new BufferedWriter(new OutputStreamWriter(f2));

            int tab_mol[] = new int[nombre];
            int mol;
            for(int i = 0; i <= nombre - 1 ; i++)
            {
                mol = Integer.parseInt(buff.readLine());
                tab_mol[i] = mol;
            }
            // compute fingerprint  of 10000 molecules
            IBitFingerprint [] tabfinger = new IBitFingerprint[nombre];
            IChemObjectBuilder builder;
            boolean skipErrors;
            String fichier1,fichier2;
            long debut,fin;
            long total1,total2;

            total1 = System.nanoTime() ;
            for(int i = 0; i < nombre; i++)
            {
                //System.out.println(i);
                try {
                    fichier1 = "CHEBI/" + tab_mol[i] + "/" + tab_mol[i] + ".sdf";
                    InputStream in = new FileInputStream(new File(fichier1));
                    Reader rdr = new InputStreamReader(in, StandardCharsets.UTF_8);
                    builder = SilentChemObjectBuilder.getInstance();
                    skipErrors = true;
                    IteratingSDFReader sdfr = new IteratingSDFReader(rdr, builder, skipErrors);

                    IAtomContainer mol1 = null;
                    while (sdfr.hasNext()) {
                        mol1 = sdfr.next();
                    }

                    debut = System.nanoTime() ;
                    //Fingerprinter fingerprint = new Fingerprinter(); the line for dAYLIGHT
                    // CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4); fingerprint ECFP4
                    CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP6);
                    tabfinger[i] = fingerprint.getBitFingerprint(mol1);
                    fin = System.nanoTime();
                    res1.write(Double.toString((fin - debut)/(double)1000000000) + "\t");
                    res1.flush();
                    res1.newLine();


                } catch (FileNotFoundException e) {
                    System.err.println("No such file: " + e.getMessage());
                    System.exit(-3);
                }


            }
            total2 = System.nanoTime() ;
            System.out.println("computation time of fingerprint: "+ (total2 - total1)/(double)1000000000);
            res1.close();

            total1 = System.nanoTime() ;

            // similarity calculation
            for(int i = 1; i <= nombre - 1; i++)
            {
                System.out.println(i);
                for(int j = 0; j < i ; j++)
                {
                    debut = System.nanoTime() ;
                    double tanimoto_coefficient = Tanimoto.calculate(tabfinger[i], tabfinger[j]);
                    fin = System.nanoTime() ;
                    res2.write(Double.toString((fin - debut)/(double)1000000000) + "\t");
                    res.write(Double.toString(tanimoto_coefficient) + "\t");
                    res.flush();

                }
                res.newLine();
                res2.newLine();
            }
            total2 = System.nanoTime() ;
            System.out.println("computation time of similarity : "+ (total2 - total1)/(double)1000000000);
            buff.close();
            res.close();
            res2.close();
            

        }
        catch (Exception e){
            System.out.println(e.toString());
            System.exit(-4);
        }
    }
}
