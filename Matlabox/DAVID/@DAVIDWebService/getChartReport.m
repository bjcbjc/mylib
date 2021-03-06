function xReturn = getChartReport(obj,args0,args1)
%getChartReport(obj,args0,args1)
%
%     Input:
%       args0 = (double)
%       args1 = (int)
%   
%     Output:
%       return{:} = ()

% Build up the argument lists.
values = { ...
   args0, ...
   args1, ...
   };
names = { ...
   'args0', ...
   'args1', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}double', ...
   '{http://www.w3.org/2001/XMLSchema}int', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://service.session.sample', ...
    'getChartReport', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getChartReport', ...
    soapMessage);
xReturn = parseSoapResponse(response);
