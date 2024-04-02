class Matrix
  attr_accessor :rows, :cols, :elements

  def initialize(rows, cols, elements)
    @rows = rows
    @cols = cols
    @elements = elements
    validate_elements
  end

  def get(i, j)
    raise ArgumentError, 'Index out of bounds' unless valid_index?(i, j)

    @elements[i][j]
  end

  def set(i, j, value)
    raise ArgumentError, 'Index out of bounds' unless valid_index?(i, j)

    @elements[i][j] = value
  end

  private

  def valid_index?(i, j)
    i.between?(0, @rows - 1) && j.between?(0, @cols - 1)
  end

  def validate_elements
    raise ArgumentError, 'Invalid number of rows' unless @elements.length == @rows

    @elements.each do |row|
      raise ArgumentError, 'Invalid number of columns' unless row.length == @cols
    end
  end
end

class Vector
  attr_accessor :dim, :elements

  def initialize(dim, elements)
    @dim = dim
    @elements = elements
    raise ArgumentError, 'Invalid number of elements' unless @elements.length == @dim
  end

  def get(i)
    raise ArgumentError, 'Index out of bounds' unless i.between?(0, @dim - 1)

    @elements[i]
  end

  def set(i, value)
    raise ArgumentError, 'Index out of bounds' unless i.between?(0, @dim - 1)

    @elements[i] = value
  end
end

class LinearAlgebra
  def self.transpose(a)
    if a.is_a?(Matrix)
      transposed_elements = a.elements.transpose
      Matrix.new(a.cols, a.rows, transposed_elements)
    elsif a.is_a?(Vector)
      Vector.new(a.dim, a.elements.clone)
    else
      raise ArgumentError, 'Invalid argument type'
    end
  end

  def self.sum(a, b)
      raise ArgumentError, 'Incompatible dimensions' unless a.rows == b.rows && a.cols == b.cols

      if a.is_a?(Matrix)
        result_elements = a.elements.map.with_index do |row, i|
          row.map.with_index { |elem, j| elem + b.elements[i][j] }
        end
        Matrix.new(a.rows, a.cols, result_elements)
      elsif a.is_a?(Vector)
        result_elements = a.elements.map.with_index { |elem, i| elem + b.elements[i] }
        Vector.new(a.dim, result_elements)
      else
        raise ArgumentError, 'Invalid argument types'
      end
    end

  def self.times(a, b)
    if a.is_a?(Numeric) && b.is_a?(Matrix)
      result_elements = b.elements.map { |row| row.map { |elem| elem * a } }
      Matrix.new(b.rows, b.cols, result_elements)
    elsif a.is_a?(Matrix) && b.is_a?(Matrix)
      raise ArgumentError, 'Incompatible dimensions' unless a.rows == b.rows && a.cols == b.cols

      result_elements = a.elements.map.with_index do |row, i|
        row.map.with_index { |elem, j| elem * b.elements[i][j] }
      end
      Matrix.new(a.rows, a.cols, result_elements)
    elsif a.is_a?(Vector) && b.is_a?(Vector)
      raise ArgumentError, 'Incompatible dimensions' unless a.dim == b.dim

      result_elements = a.elements.map.with_index { |elem, i| elem * b.elements[i] }
      Vector.new(a.dim, result_elements)
    else
      raise ArgumentError, 'Invalid argument types'
    end
  end

  def self.dot(a, b)
    raise ArgumentError, 'Incompatible dimensions' unless a.cols == b.rows

    result_elements = Array.new(a.rows) { Array.new(b.cols, 0) }

    a.rows.times do |i|
      b.cols.times do |j|
        a.cols.times do |k|
          result_elements[i][j] += a.get(i, k) * b.get(k, j)
        end
      end
    end

    Matrix.new(a.rows, b.cols, result_elements)
  end

  def self.gauss(a)
    raise ArgumentError, 'Matrix must be square' unless a.rows == a.cols

    n = a.rows
    m = a.elements.map(&:clone)

    (0...n).each do |i|
      pivot = m[i][i].to_f
      raise ArgumentError, 'Matrix is singular' if pivot.zero?

      (i + 1...n).each do |k|
        factor = m[k][i] / pivot
        (i...n).each { |j| m[k][j] -= factor * m[i][j] }
      end
    end

    Matrix.new(n, n, m)
  end

  def self.solve(a)
    gauss_result = gauss(a)
    n = gauss_result.rows
    x = Array.new(n, 0)

    (n - 1).downto(0) do |i|
      x[i] = gauss_result.get(i, n) / gauss_result.get(i, i)
      (i - 1).downto(0) do |k|
        gauss_result.set(k, n, gauss_result.get(k, n) - gauss_result.get(k, i) * x[i])
      end
    end

    Vector.new(n, x)
  end
end

class LinearAlgebraInterface
  def initialize
    puts "Bem-vindo à Interface de Álgebra Linear!"
    puts "Selecione uma opção:"
    puts "1. Transpor uma matriz ou vetor"
    puts "2. Somar duas matrizes ou vetores"
    puts "3. Multiplicar duas matrizes ou vetores"
    puts "4. Multiplicar uma matriz por um escalar"
    puts "5. Realizar o produto escalar de dois vetores"
    puts "6. Resolver um sistema de equações lineares"
    puts "7. Realizar a eliminação gaussiana"
    puts "8. Sair"

    option = gets.chomp.to_i

    case option
    when 1
      transpose_interface
    when 2
      sum_interface
    when 3
      times_interface
    when 4
      scalar_multiplication_interface
    when 5
      dot_product_interface
    when 6
      solve_equation_interface
    when 7
      gaussian_elimination_interface
    when 8
      exit
    else
      puts "Opção inválida. Por favor, selecione uma opção válida."
      initialize
    end
  end

  private

  def times_interface
      puts "Multiplicar duas matrizes ou vetores"

      puts "Digite 1 para multiplicar duas matrizes, 2 para multiplicar uma matriz por um vetor ou 3 para multiplicar dois vetores:"
      type = gets.chomp.to_i

      case type
      when 1
        puts "Digite a primeira matriz:"
        matrix1 = get_matrix_input
        puts "Digite a segunda matriz:"
        matrix2 = get_matrix_input
        result = LinearAlgebra.times(matrix1, matrix2)
        puts "Resultado da multiplicação:"
        print_matrix(result)
      when 2
        puts "Digite a matriz:"
        matrix = get_matrix_input
        puts "Digite o vetor:"
        vector = get_vector_input
        result = LinearAlgebra.times(matrix, vector)
        puts "Resultado da multiplicação:"
        print_vector(result)
      when 3
        puts "Digite o primeiro vetor:"
        vector1 = get_vector_input
        puts "Digite o segundo vetor:"
        vector2 = get_vector_input
        result = LinearAlgebra.times(vector1, vector2)
        puts "Resultado da multiplicação:"
        print_vector(result)
      else
        puts "Opção inválida. Por favor, selecione 1, 2 ou 3."
        times_interface
      end
    end

  def get_matrix_input
      print "Digite o número de linhas: "
      rows = gets.chomp.to_i
      print "Digite o número de colunas: "
      cols = gets.chomp.to_i
      elements = []

      rows.times do |i|
        print "Digite os elementos da linha #{i + 1} separados por espaço: "
        row_elements = gets.chomp.split.map(&:to_f)
        elements << row_elements
      end

      Matrix.new(rows, cols, elements)
    end

    def get_vector_input
      print "Digite a dimensão do vetor: "
      dim = gets.chomp.to_i
      print "Digite os elementos do vetor separados por espaço: "
      elements = gets.chomp.split.map(&:to_f)

      Vector.new(dim, elements)
    end


    def transpose_interface
      puts "Transpor uma matriz ou vetor"

      puts "Digite 1 para uma matriz e 2 para um vetor:"
      type = gets.chomp.to_i

      case type
      when 1
        matrix = get_matrix_input
        result = LinearAlgebra.transpose(matrix)
        puts "Matriz transposta:"
        print_matrix(result)
      when 2
        vector = get_vector_input
        result = LinearAlgebra.transpose(vector)
        puts "Vetor transposto:"
        print_vector(result)
      else
        puts "Opção inválida. Por favor, selecione 1 ou 2."
        transpose_interface
      end
    end

    def sum_interface
      puts "Somar duas matrizes ou vetores"

      puts "Digite 1 para duas matrizes e 2 para dois vetores:"
      type = gets.chomp.to_i

      case type
      when 1
        puts "Digite a primeira matriz:"
        matrix1 = get_matrix_input
        puts "Digite a segunda matriz:"
        matrix2 = get_matrix_input
        result = LinearAlgebra.sum(matrix1, matrix2)
        puts "Resultado da soma:"
        print_matrix(result)
      when 2
        puts "Digite o primeiro vetor:"
        vector1 = get_vector_input
        puts "Digite o segundo vetor:"
        vector2 = get_vector_input
        result = LinearAlgebra.sum(vector1, vector2)
        puts "Resultado da soma:"
        print_vector(result)
      else
        puts "Opção inválida. Por favor, selecione 1 ou 2."
        sum_interface
      end
    end

    def dot_product_interface
      puts "Realizar o produto escalar de dois vetores"

      puts "Digite o primeiro vetor:"
      vector1 = get_vector_input

      puts "Digite o segundo vetor:"
      vector2 = get_vector_input

      result = LinearAlgebra.dot(vector1, vector2)

      puts "Resultado do produto escalar:"
      puts result.elements
    end

    def solve_equation_interface
      puts "Resolver um sistema de equações lineares"

      puts "Digite a matriz de coeficientes:"
      coefficients_matrix = get_matrix_input

      puts "Digite o vetor de termos constantes:"
      constants_vector = get_vector_input

      result = LinearAlgebra.solve(Matrix.new(coefficients_matrix.rows, coefficients_matrix.cols + 1, coefficients_matrix.elements.map.with_index { |row, i| row + [constants_vector.elements[i]] }))

      puts "Solução do sistema de equações:"
      puts result.elements
    end

    def scalar_multiplication_interface
      puts "Multiplicar uma matriz por um escalar"

      puts "Digite o escalar:"
      scalar = gets.chomp.to_f

      puts "Digite a matriz:"
      matrix = get_matrix_input

      result = LinearAlgebra.times(scalar, matrix)

      puts "Resultado da multiplicação:"
      print_matrix(result)
    end



    def gaussian_elimination_interface
      puts "Realizar a eliminação gaussiana"

      puts "Digite a matriz:"
      matrix = get_matrix_input

      result = LinearAlgebra.gauss(matrix)

      puts "Matriz após a eliminação gaussiana:"
      print_matrix(result)
    end


    # Métodos similares para as outras operações...

    def print_matrix(matrix)
      matrix.elements.each do |row|
        puts row.join(" ")
      end
    end

    def print_vector(vector)
      puts vector.elements.join(" ")
    end
  end

  LinearAlgebraInterface.new
