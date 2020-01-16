using System;

namespace DiplomaProject
{
    public static class Display
    {
        public static void Array<T>(T[] array)
        {
            foreach (var element in array)
            {
                Console.Write($"{element} ");
            }
            Console.WriteLine();
        }

        public static void Matrix<T>(T[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    Console.Write($"{matrix[i,j]:0.####} \t");
                }
                Console.WriteLine();
            }
        }
    }
}