package MyMath;

/**
 * Created by nikitafedorov on 29/03/2017.
 */
public class Factorial {
    public static long get(int value){
        long answer = 1;
        for(int i = 1; i <= value; i++){
            answer *= i;
        }
        return answer;
    }
}
