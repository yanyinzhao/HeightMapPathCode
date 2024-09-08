void hash_function_two_keys_to_one_key(unsigned x_vertex_num, unsigned x, unsigned y, unsigned &x_y)
{
    x_y = y * x_vertex_num + x;
}

void hash_function_two_keys_to_one_key_int(int row_or_column, int i, int j, int &i_j)
{
    i_j = i * row_or_column + j;
}

void hash_function_one_key_to_two_keys(unsigned x_vertex_num, unsigned x_y, unsigned &x, unsigned &y)
{
    x = x_y % x_vertex_num;
    y = x_y / x_vertex_num;
}