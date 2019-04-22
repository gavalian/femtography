import org.jlab.jnp.readers.*;
import org.jlab.groot.data.*;
import org.jlab.groot.ui.*;

String filename = args[0];


String output_R = filename + "_REAL.nrrd";
String output_I = filename + "_IM.nrrd";

TextFileReader reader = new TextFileReader();
reader.open(filename);

int xBin = 0;
int tBin = 0;

H2F Hr = new H2F("Hr",200,0,1.0,200,0.0,1.0);
H2F Hi = new H2F("Hi",200,0,1.0,200,0.0,1.0);
TCanvas c = new TCanvas("c",900,450);

while(reader.readNext()==true){
  double value_r = reader.getAsDouble(2);
  double value_i = reader.getAsDouble(4);
  System.out.println(tBin + " : " + xBin + " = " + value_r );
  Hr.setBinContent(tBin,xBin,value_r);
  Hi.setBinContent(tBin,xBin,value_i);


  //if(xBin>20) H.setBinContent(tBin,xBin-20,value);
  tBin++;
  if(tBin==200){
    tBin = 0; xBin++;
  }
}

H3F H3_R = Hr.createH3F(200);
H3F H3_I = Hi.createH3F(200);

H3_R.export(output_R);
H3_I.export(output_I);

c.divide(2,1);
c.cd(0).draw(Hr);
c.cd(1).draw(Hi);
