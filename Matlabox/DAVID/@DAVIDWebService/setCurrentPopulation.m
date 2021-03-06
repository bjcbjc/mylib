function xReturn = setCurrentPopulation(obj,args0)
%setCurrentPopulation(obj,args0)
%
%     Input:
%       args0 = (int)
%   
%     Output:
%       return = (int)

% Build up the argument lists.
values = { ...
   args0, ...
   };
names = { ...
   'args0', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}int', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://service.session.sample', ...
    'setCurrentPopulation', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:setCurrentPopulation', ...
    soapMessage);
xReturn = parseSoapResponse(response);
