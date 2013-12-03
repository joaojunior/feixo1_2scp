execute PARAMS {
  cplex.tilim = 7200;
}

int number_items = ...;
int costs[0..number_items-1] = ...;
int weights[0..number_items-1] = ...;  
int weight_max = ...;
dvar int+ x[0..number_items-1];

minimize
    sum(i in 0..number_items-1)
		costs[i] * x[i];
subject to{
  	sum(i in 0..number_items-1) weights[i] * x[i] <= weight_max;
};