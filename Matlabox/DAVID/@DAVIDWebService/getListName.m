function xReturn = getListName(obj,args0,args1)
%getListName(obj,args0,args1)
%
%     Input:
%       args0 = (int)
%       args1 = (int)
%   
%     Output:
%       return = (string)

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
   '{http://www.w3.org/2001/XMLSchema}int', ...
   '{http://www.w3.org/2001/XMLSchema}int', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://service.session.sample', ...
    'getListName', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getListName', ...
    soapMessage);
xReturn = parseSoapResponse(response);
