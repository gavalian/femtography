import org.jlab.jnp.readers.*;
import org.jlab.groot.data.*;
import org.jlab.groot.ui.*;

String filename = args[0];

TextFileReader reader = new TextFileReader();
reader.open(filename);

int xBin = 0;
int tBin = 0;

H2F H = new H2F("H",200,0,1.0,180,0.0,1.0);
TCanvas c = new TCanvas("c",500,500);

while(reader.readNext()==true){
  double value = reader.getAsDouble(2);
  System.out.println(tBin + " : " + xBin + " = " + value );
  if(xBin>20) H.setBinContent(tBin,xBin-20,value);
  tBin++;
  if(tBin==200){
    tBin = 0; xBin++;
  }
}

double hmin = H.getMin();
double hmax = H.getMax();
System.out.println("min/max = " + hmin + " / " + hmax);

c.draw(H);
