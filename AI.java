import java.util.function.Function;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.Random;
import java.util.Set;

/**
 * @author Ben Russell Martin
 * @version 09/11/2000
 * @see GA
 */
public class AI {

    /**
     * 
     * @param args
     */
    public static void main(String[] args){
        //Do not delete/alter the next line
        long startT=System.currentTimeMillis();

        // Create a the object for problem one with a population size of 50, a comparator favoring lower fitnesses, and a supplied fitness function
        Problem1 problem1 = new Problem1(50, (x1, x2) -> Double.valueOf(x1.getFitness()).compareTo(x2.getFitness()), Assess::getTest1);
        Thread t_problem1 = new Thread(problem1); // Create a thread for the problem
        
        // Create a the object for problem one with a population size of 150, a comparator favoring higher fitnesses, and a supplied fitness function
        Problem2 problem2 = new Problem2(150, (x1, x2) -> Double.valueOf(x2.getFitness()).compareTo(x1.getFitness()), Assess::getTest2);
        Thread t_problem2 = new Thread(problem2); // Create a thread for the problem

        // Start both solutions implementing parallelisation
        t_problem1.start();
        t_problem2.start();

        try {
            t_problem1.join();
            t_problem2.join();
            // Wait until both solutions have completed
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        // Fetch and check-in the proposed candidate solutions for assessment
        double[] solution1 = problem1.getBest().getVals();
        boolean[] solution2 = problem2.getBest().getVals();
        // Submit solutions
        Assess.checkIn("Ben Martin", "bm433", solution1, solution2);

        // Display execution time for stats
        System.out.printf("T: %f%n", ((System.currentTimeMillis() - startT) / 1000.0));
    }
    
    /**
     * A base class used to generate the "Problem1" and "Probelm2" classes
     */
    private abstract static class GA<T, F> implements Runnable {
        // T = candidate solution encoding type
        // F = fitness function return type
        protected final int POPULATIONSIZE;
        protected static Random r = new Random();
        protected Set<CandidateSolution> pop;
        protected Function<T, F> f;
        private Comparator<CandidateSolution> c;
        
        /**
         * Create a genetic algorighm object with the population size, the fitness comparator, and fitness funciton handed as params
         * @param populationSize The fixed population size of each generation of the type {@link int}.
         * @param c The comparator used to check which candidate solution has the better fitness of the type {@link Comparator}
         * @param f The fitness function used to generate the fitness coefficient for a particular {@link CandidateSolution} of the type {@link Function}.
         */
        public GA(int populationSize, Comparator<CandidateSolution> c, Function<T, F> f) {
            // Create the state
            this.POPULATIONSIZE = populationSize;
            this.pop = new HashSet<>();
            this.c = c;
            this.f = f;
        }

        /**
         * Return a {@link String} containing a list of fitnesses for the current generation population.
         * Usefull for debbugging
         */
        @Override
        public String toString() {
            return Arrays.toString(sortPop().stream().mapToDouble(CandidateSolution::getFitness).toArray());
        }

        /**
         * Return a sorted list of {@link CandidateSolutions} by their fitness coefficient using the {@link Comparator} in the {@code c} field.
         * Sorted by best to worse.
         * @return Sorted list of {@link CandidateSolution}s as a {@link List}.
         */
        public List<CandidateSolution> sortPop() {
            return pop.stream().sorted(c).collect(Collectors.toList());
        }

        /**
         * Return the best {@link CandidateSolution} from the current generation population.
         * @return The best candidate solution as a {@link CandidateSolution}
         */
        public CandidateSolution getBest() {
            return sortPop().get(0);
        }
        
        /**
         * Method containing the generation loop and the genetic algoithm control flow
         * @return The encoding of the best candidate solution
         */
        public abstract void run();

        /**
         * A base class used to generate the candidate solutions.
         */
        public abstract class CandidateSolution {
            protected T vals;

            /**
             * Construct a candidate solution with the encoding passed to its params.
             * @param vals The encoding of the new candidate solution
             */
            protected CandidateSolution(T vals) {
                this.vals = vals;
            }

            /**
             * Getter for the encoding of the candidate solution
             * @return The encoding of the candidate solution.
             */
            public T getVals() {
                return vals;
            }

            /**
             * Getter method for the fitness
             * @return The fitness of the candidate solution as a {@code double}.
             */
            public abstract double getFitness();
            
            /**
             * Implements the method of mutating the candidate solution
             * @param maxFactor The factor by which the candidate solution gets mutated
             * @return The new mutated candidate solution as a {@link CandidateSolution}
             */
            public abstract CandidateSolution mutate(double maxFactor);

            /**
             * Implements the method of crossover for two {@link CandidateSolution} parents - itself and one other handed as a param.
             * @param other The other parent for the crossover as a {@link CandidateSolution}.
             * @return The new mutated candidate solution child as a {@link CandidateSolution}
             */
            public abstract CandidateSolution crossover(CandidateSolution other);
        }

    }
    
    public static class Problem1 extends GA<double[], Double> {
        private double truncationRatio = .1; // Ratio used in the truncation selection method https://en.wikipedia.org/wiki/Truncation_selection

        public Problem1(int populationSize, Comparator<CandidateSolution> c, Function<double[], Double> f) {
            super(populationSize, c, f);
        }

        /**
         * Helper method for generating a new population
         * @return a {@link Set} of {@link CandidateSolution}s containing a new populaton.
         */
        protected Set<CandidateSolution> generatePop() {
            Set<CandidateSolution> newPop = new HashSet<>();
            for (int i = 0; i < POPULATIONSIZE; i++) {
                // Iterate until the population size requirement is met generating a new candidate solution per iteration
                // Generate a new anonymous object of the type CandidateSolution and add it to the the new population
                newPop.add(new CandidateSolution(Arrays.stream(new double[20]).map(x -> (r.nextDouble() * 10) - 5).toArray()) {
                    public CandidateSolution crossover(CandidateSolution other) {
                        double[] otherVals = other.getVals();
                        double[] newVals = new double[vals.length];

                        // Splice
                        // Generate a random integer taking all encodings from the first parent up (itslef) to the random integer and the rest of the encodings from the other parent (handed to as a method param).
                        int spliceIndex = r.nextInt(vals.length);
                        for (int i = 0; i < spliceIndex; i++) {
                            newVals[i] = vals[i];
                        }
                        for (int i = spliceIndex; i < vals.length; i++) {
                            newVals[i] = otherVals[i];
                        }

                        try {
                            // https://docs.oracle.com/javase/specs/jls/se15/html/jls-15.html#jls-15.9.3
                            // Generate a new instance of the anonymous object with the new encoding and return
                            return getClass().getDeclaredConstructor(Problem1.class, double[].class).newInstance(Problem1.this, newVals);
                        } catch (Exception e) {
                            // Generic catch for the getDeclaredConstructor and newInstance throws.
                            e.printStackTrace();
                            return null;
                        }
                    }

                    public CandidateSolution mutate(double maxFactor) {
                        double[] newVals = vals;
                        // Add to a selected element in the encoding by a random amount between +maxFactor and -maxFactor
                        int mutateIndex = r.nextInt(newVals.length);
                        newVals[mutateIndex] += (r.nextDouble() * 2 * maxFactor) - maxFactor;
                        try {
                            // https://docs.oracle.com/javase/specs/jls/se15/html/jls-15.html#jls-15.9.3
                            // Generate a new instance of the anonymous object with the new encoding and return
                            return getClass().getDeclaredConstructor(Problem1.class, double[].class).newInstance(Problem1.this, newVals);
                        } catch (Exception e) {
                            // Generic catch for the getDeclaredConstructor and newInstance throws.
                            e.printStackTrace();
                            return null;
                        }
                    }

                    public double getFitness() {
                        // Fetch the fitness by calling the fitness function and returning its result.
                        return f.apply(vals);
                    }

                });
            }
            // Return the newly generated population
            return newPop;
        }

        @Override
        public void run() {
            int generations = 0; // Generation counter for performance statistics
            long t0 = System.nanoTime(); // Starting time for performance statistcs
            pop = generatePop(); // Starting population
            int truncationPop = (int) Math.floor(POPULATIONSIZE * truncationRatio); // The size of the truncated population
            double bestFitness;
            do {
                List<CandidateSolution> newPop = sortPop().subList(0, truncationPop); // Population for the next generation
                CandidateSolution best = getBest();
                bestFitness = best.getFitness();
                while (newPop.size() < POPULATIONSIZE) {
                    // Crossover and mutate from the truncated population until the new population meets population size requirements
                    CandidateSolution candidateSolution = newPop.get(r.nextInt(truncationPop)).crossover(newPop.get(r.nextInt(truncationPop)));
                    candidateSolution = candidateSolution.mutate((bestFitness * 1) + .0000001);
                    newPop.add(candidateSolution);
                }
                // Replace the current population with the new population
                pop = new HashSet<>(newPop);
                generations++;

            } while (System.nanoTime() - t0 < 5e+9 && bestFitness != 0); // Iterate until either the GA takes too long or the the fitness is perfect
            
            // Display statistics
            System.out.printf("Problem 1: %e fitness \u2014 %d cycles in %fs%n", getBest().getFitness(), generations, (System.nanoTime() - t0) / 1e+9);
        }
        
    }

    public static class Problem2 extends GA<boolean[], double[]> {

        public Problem2(int populationSize, Comparator<CandidateSolution> c, Function<boolean[], double[]> f) {
            super(populationSize, c, f);
        }

        /**
         * Helper method for generating a new population
         * @return a {@link Set} of {@link CandidateSolution}s containing a new populaton.
         */
        protected Set<CandidateSolution> generatePop() {
            // Create set to hold the population and fill it with candidate solution objects
            Set<CandidateSolution> newPop = new HashSet<>();
            for (int i = 0; i < POPULATIONSIZE; i++) {
                boolean[] vals = new boolean[100];
                for (int j = 0; j < vals.length; j++) {
                    vals[j] = false;
                    // Start with all candidate solutions containing nothing in the bag
                    // Prevevents the chances of all candidate solutions having a weight over the max forfeiting a generation
                }

                newPop.add(new CandidateSolution(vals) {
                    public CandidateSolution crossover(CandidateSolution other) {
                        boolean[] otherVals = other.getVals();
                        boolean[] newVals = new boolean[vals.length];

                        // Splice
                        // Similar to problem 1
                        // Chose an index to splice at random and joining the left elements of one parent with the right of another.
                        int spliceIndex = r.nextInt(vals.length);
                        for (int i = 0; i < spliceIndex; i++) {
                            newVals[i] = vals[i];
                        }
                        for (int i = spliceIndex; i < vals.length; i++) {
                            newVals[i] = otherVals[i];
                        }

                        try {
                            // https://docs.oracle.com/javase/specs/jls/se15/html/jls-15.html#jls-15.9.3
                            // Create a new instance of the candidate solution with the new values
                            return getClass().getDeclaredConstructor(Problem2.class, boolean[].class).newInstance(Problem2.this, newVals);
                        } catch (Exception e) {
                            e.printStackTrace();
                            return null;
                        }
                    }

                    public CandidateSolution mutate(double maxFactor) {
                        // Bit string mutation
                        boolean[] newVals = vals;
                        for (int j = 0; j < vals.length * maxFactor; j++) {
                            // Flip the bits of maxFactor% of the encoding
                            // i.e: 1 becomes 0 and 0 becomes 1 of a randomly chosen bits in the encoding
                            // More than one chosen to prevent getting stuck in local minimas
                            int mutateIndex = r.nextInt(newVals.length);
                            newVals[mutateIndex] = !newVals[mutateIndex];
                        }
                        try {
                            // Create a new instance of the candidate solution with the new values
                            return getClass().getDeclaredConstructor(Problem2.class, boolean[].class).newInstance(Problem2.this, newVals);
                        } catch (Exception e) {
                            e.printStackTrace();
                            return null;
                        }
                    }

                    public double getFitness() {
                        // Return 0 if the weight is greater than the max, otherwise return the weight
                        double[] fitness = f.apply(vals);
                        if (fitness[0] > 500) {
                            return 0;
                        }
                        return fitness[1];
                    }

                });
            }
            return newPop;
        }

        @Override
        public void run() {
            // For statistics
            long t0 = System.nanoTime();
            int generations = 0;

            pop = generatePop(); // Create the starting generation
            double truncationRatio = .1; // The ratio of the first generation to bring to the second during selection
            int truncationPop = (int) Math.floor(POPULATIONSIZE * truncationRatio); // The count of the number of candidates to bring to the second generation
            do {
                List<CandidateSolution> newPop = sortPop().subList(0, truncationPop); // Get the truncated population for the next generation
                for (int i = POPULATIONSIZE - truncationPop; i < POPULATIONSIZE; i++) {
                    // Create new candidate solutions until population size requirement met
                    newPop.add(
                        newPop.get(r.nextInt(truncationPop)) // Get a member of the previous generation from the truncated population
                        .crossover(newPop.get(r.nextInt(truncationPop))) // Crossover the selected candidate solution with another member of the truncated population
                        .mutate(0.02) // Mutate the crossedover candidate solution by a factor of 2%
                    );
                }
                // Replace the popuation with the new candidate solutions and increment generation counter for stats
                pop = new HashSet<>(newPop);
                generations++;
            } while (System.nanoTime() - t0 < 5e+9 && getBest().getFitness() < 203); // Run until out of time or fitness is above 29 of the performance marks (fitness > 203)
            
            // Display statistics
            System.out.printf("Problem 2: %f fitness \u2014 %d cycles in %fs%n", getBest().getFitness(), generations, (System.nanoTime() - t0) / 1e+9);
        }

    }
    
}
