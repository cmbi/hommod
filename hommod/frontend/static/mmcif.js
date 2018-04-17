function get_cif_values(s)
{
  s = s.trim();
  var values = [];
  while (s.length > 0)
  {
    if (s.charAt(0) == '\'' || s.charAt(0) == '\"')
    {
      var quote = s.charAt(0);
      var i = s.substring(1).search(quote)
      if (i == -1)
        throw new SyntaxError("Closing " + quote + " expected in " + s);
      values.push(s.substring(1, i + 1));
      s = s.substring(i + 2).trim();
    }
    else
    {
      var v = s.split(/\s+/)[0];
      values.push(v);
      s = s.substring(v.length).trim();
    }
  }
  return values;
}

function parse_mmcif(content)
{
  var tables = {};

  var current_row = {},
      current_table_id = '',
      values = [],
      header = [],
      inloop = false;

  var lines = content.split("\n");
  for (var i = 0; i < lines.length; i++)
  {
    var line = lines[i].trim();

    if (line.startsWith('data_') || line.startsWith('#'))
    {
      if (values.length > 0 && values.length < header.length)
        throw new SyntaxError("Too few values in " + values + "(" + values.length + ")" + ", " + header.length + " expected");

      if (Object.keys(current_row).length > 0)
        tables[current_table_id].push(current_row);

      inloop = false;
      header = [];
      values = [];
      current_row = {};
    }
    else if (line == "loop_")
      inloop = true;
    else if (line.startsWith('_'))  // Starts with a variable id
    {
      var id = line.split(/\s+/)[0];
      var s = id.split('.');
      current_table_id = s[0];
      var var_id = s[1];

      if (!(current_table_id in tables))
        tables[current_table_id] = [];

      if (inloop)
        header.push(var_id);
      else
      {
        var v = line.substring(id.length).trim();
        if (v.length == 0)
        {
          i++;
          line = lines[i].trim();

          if (line.startsWith(';'))  // multiline value
          {
            v = line.substring(1).trim();

            while (true)
            {
              i++;
              line = lines[i].trim();
              if (line.startsWith(';'))
                break;
              v += line.trim();
            }
            values = [v];
          }
          else
            values = get_cif_values(line);
        }
        else
          values = get_cif_values(v);

        if (values.length > 1)
          throw new SyntaxError("only 1 value expected for " + current_table_id + "." + var_id + ", got: " + values);

        current_row[var_id] = values[0];
        values = [];
      }
    }
    else if (inloop)  // No variable id expected on this line
    {
      if (line.startsWith(';'))  // multiline value
      {
        var j = values.length;
        values.push(line.substring(1).trim());

        while (true)
        {
          i++;
          line = lines[i].trim();
          if (line.startsWith(';'))
            break;
          values[j] += line;
        }
      }
      else
        values = values.concat(get_cif_values(line));

      if (values.length == header.length)
      {
        current_row = {};
        for (var j = 0; j < header.length; j++)
          current_row[header[j]] = values[j];
        tables[current_table_id].push(current_row);
        current_row = {};
        values = [];
      }
      else if (values.length > header.length)
        throw new SyntaxError("Too many values in " + values + "(" + values.length + ")" + ", " + header.length + " expected");
    }
  }

  return tables;
}
