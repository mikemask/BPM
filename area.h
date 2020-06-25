 double area(double dist, double &d1, double &d2)
{ 
  double R1 = d1/2.;
  double R2 = d2/2.;
  double a =3.14159265*pow((R1*sin(acos((pow(R1,2.)-pow(R2,2)+pow(dist,2.)/(2*R1*dist))))),2.);
  return a;
} 

